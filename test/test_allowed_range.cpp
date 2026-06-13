// allowed_range の結果（tquad 内側の hmesh エッジ集合）を期待値 JSON と突き合わせる回帰テスト。
// polyscope 不要。期待値は別ファイル expected_allowed.json に持つ（再ビルド不要で編集可）。
//
// JSON 形式（このテストが必要とするサブセット）:
//   {
//     "<mesh basename>": { "<tqid>": [eid, eid, ...], ... },
//     ...
//   }
//   ・配列 = その tquad で内側であるべき hmesh エッジ id（順不同。重複は無視）。
//   ・配列が空 [] の tquad は「期待値未設定」とみなしスキップし、計算結果を [fill] 行として出力する
//     （可視化で正しいと確認できたら、その出力を JSON に貼り付けて確定させる用途）。
//
// 依存を増やさないため、上記サブセット専用の小さな JSON パーサを内蔵する。
#include <igl/readOBJ.h>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cctype>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "check.h"

using namespace metriko;
using Expected = std::map<std::string, std::map<std::string, std::vector<int>>>;

// ---- 最小 JSON サブセットパーサ: { "k": { "k": [int,...] } } のみ対応 ----
struct JParser {
    const std::string& s;
    size_t i = 0;
    bool ok = true;

    void ws() { while (i < s.size() && std::isspace((unsigned char)s[i])) ++i; }
    bool eat(char c) { ws(); if (i < s.size() && s[i] == c) { ++i; return true; } return false; }
    std::string str() {
        ws();
        if (i >= s.size() || s[i] != '"') { ok = false; return {}; }
        ++i; std::string r;
        while (i < s.size() && s[i] != '"') r += s[i++];
        if (i >= s.size()) { ok = false; return {}; }
        ++i; return r;
    }
    int integer() {
        ws(); size_t j = i;
        while (i < s.size() && (std::isdigit((unsigned char)s[i]) || s[i] == '-' || s[i] == '+')) ++i;
        if (i == j) { ok = false; return 0; }
        return std::stoi(s.substr(j, i - j));
    }
    std::vector<int> array() {
        std::vector<int> v;
        if (!eat('[')) { ok = false; return v; }
        if (eat(']')) return v;
        do { v.push_back(integer()); } while (eat(','));
        if (!eat(']')) ok = false;
        return v;
    }
    // { "tqid": [..], ... }
    std::map<std::string, std::vector<int>> inner() {
        std::map<std::string, std::vector<int>> m;
        if (!eat('{')) { ok = false; return m; }
        if (eat('}')) return m;
        do { std::string k = str(); if (!eat(':')) { ok = false; break; } m[k] = array(); } while (eat(','));
        if (!eat('}')) ok = false;
        return m;
    }
    Expected top() {
        Expected e;
        if (!eat('{')) { ok = false; return e; }
        if (eat('}')) return e;
        do { std::string k = str(); if (!eat(':')) { ok = false; break; } e[k] = inner(); } while (eat(','));
        if (!eat('}')) ok = false;
        return e;
    }
};

static std::string basename_we(const std::string& p) {
    auto s = p.find_last_of("/\\");
    std::string f = (s == std::string::npos) ? p : p.substr(s + 1);
    auto d = f.find_last_of('.');
    return (d == std::string::npos) ? f : f.substr(0, d);
}

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: test_allowed_range <gridscale> <expected.json|--dump> <mesh.obj> [more.obj ...]\n"; return 2; }
    const int N = 4;
    const double scale = std::stod(argv[1]);

    // --dump: 全 tquad の allowed_range_trace 結果を JSON で出力（golden 採取用）
    const bool dump = std::string(argv[2]) == "--dump";

    Expected expected;
    if (!dump) {
        std::ifstream jf(argv[2]);
        CHECK(jf.good());                              // JSON が開けること
        std::stringstream ss; ss << jf.rdbuf();
        std::string text = ss.str();
        JParser jp{ text };
        expected = jp.top();
        CHECK(jp.ok);                                  // JSON が（サブセットとして）パースできること
    }

    int n_checked = 0, n_skipped = 0;
    if (dump) std::cout << "{\n";

    for (int a = 3; a < argc; ++a) {
        const char* mesh = argv[a];
        const std::string key = basename_we(mesh);
        if (!dump && !expected.count(key)) { std::cout << "[skip] no expected entry for '" << key << "'\n"; continue; }

        // --- pipeline: field -> param -> mc -> Tmesh -> TmeshMut（test_tmesh と同一）---
        MatXd V; MatXi F;
        igl::readOBJ(mesh, V, F);
        CHECK(V.rows() > 0 && F.rows() > 0);
        Hmesh hm(V, F);

        FaceRosyField rawf(hm, N, FieldType::Smoothest);
        rawf.computeMatching(MatchingType::Principal);
        auto seam = compute_seam(rawf);
        auto cutm = compute_cut_mesh(hm, seam);
        auto cmbf = compute_combbed_field(rawf, seam);

        MatXd ext(hm.nF, 3 * N);
        for (Face f : hm.faces)
            for (int k = 0; k < N; ++k) {
                complex c = cmbf->field(f.id, k);
                ext.block(f.id, 3 * k, 1, 3) = (c.real() * f.basisX() + c.imag() * f.basisY()).normalized();
            }

        RosyParameterization rp(hm, *cutm, ext, cmbf->singular, cmbf->matching, seam, N, scale);
        rp.seamless = false;
        rp.localInjectivity = true;
        rp.verbose = false;
        rp.setup();
        rp.integ();

        VecXc uv2(hm.nF * 3);
        for (Face f : hm.faces) {
            uv2(f.id * 3 + 0) = complex{ rp.cfn(f.id, 0), rp.cfn(f.id, 1) };
            uv2(f.id * 3 + 1) = complex{ rp.cfn(f.id, 4), rp.cfn(f.id, 5) };
            uv2(f.id * 3 + 2) = complex{ rp.cfn(f.id, 8), rp.cfn(f.id, 9) };
        }

        mc::Mgrph mg(hm, uv2, cmbf->matching, cmbf->singular);
        Tmesh    tm(mg);
        TmeshMut tmm(mg, tm);

        if (dump) {                                    // 全 tquad を JSON 出力して次のメッシュへ
            std::cout << "  \"" << key << "\": {";
            for (int tqid = 0; tqid < (int)tmm.tquads.size(); ++tqid) {
                std::set<int> got;
                for (auto& [eid, r0, r1] : tmm.allowed_range_trace(tqid)) got.insert(eid);
                std::cout << (tqid ? "," : "") << "\n    \"" << tqid << "\": [";
                bool first = true;
                for (int e : got) { std::cout << (first ? "" : ",") << e; first = false; }
                std::cout << "]";
            }
            std::cout << "\n  }" << (a + 1 < argc ? "," : "") << "\n";
            continue;
        }

        // --- 各 tqid について allowed_range を期待値と比較 ---
        for (auto& [tqid_str, exp_vec] : expected[key]) {
            int tqid = std::stoi(tqid_str);
            CHECK(tqid >= 0 && tqid < (int)tmm.tquads.size());

            std::set<int> got;
            for (auto& [eid, r0, r1] : tmm.allowed_range_trace(tqid)) got.insert(eid);

            std::set<int> exp(exp_vec.begin(), exp_vec.end());

            if (exp.empty()) {                         // 期待値未設定 → 計算結果を出力してスキップ
                std::cout << "[fill] \"" << key << "\": { \"" << tqid << "\": [";
                bool first = true;
                for (int e : got) { std::cout << (first ? "" : ",") << e; first = false; }
                std::cout << "] }\n";
                ++n_skipped;
                continue;
            }

            if (got != exp) {
                std::cerr << "FAIL: allowed_range mismatch  " << key << " tq" << tqid
                          << "  expected=" << exp.size() << " got=" << got.size() << "\n";
                auto dump = [](const char* tag, const std::set<int>& a, const std::set<int>& b) {
                    std::cerr << "  " << tag << ":";
                    int n = 0;
                    for (int e : a) if (!b.count(e)) { std::cerr << " " << e; if (++n >= 40) { std::cerr << " ..."; break; } }
                    std::cerr << "\n";
                };
                dump("missing (expected but not produced)", exp, got);
                dump("extra   (produced but not expected)", got, exp);
                return 1;
            }
            std::cout << "[ok]   " << key << " tq" << tqid << "  edges=" << got.size() << "\n";
            ++n_checked;
        }
    }

    if (dump) { std::cout << "}\n"; return 0; }

    std::cout << "[test_allowed_range] checked=" << n_checked << " skipped(unset)=" << n_skipped << "\n";
    return 0;
}
