// Regression check: compare allowed_range(tqid) of every tquad against a saved
// fixture snapshot to detect degradation.
//
// fixture format (test/fixtures/tquad_region/<mesh>.json, produced by dump_region):
//   { "<tqid>": [[eid, r0, r1], [eid, r0, r1], ...], ... }
//
// usage: test_allowed_range <gridscale> <mesh.obj> <fixture.json>
//   exits non-zero if any tquad's allowed_range differs from the fixture.
#include <string>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include "pipeline.h"
#include "check.h"

using namespace metriko;

using Tri  = std::tuple<int, double, double>;       // (eid, r0, r1)
using Tris = std::vector<Tri>;

// Minimal JSON parser for { "<tqid>": [[int,double,double], ...] }.
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
    double num() {
        ws(); size_t j = i;
        while (i < s.size() && (std::isdigit((unsigned char)s[i]) ||
                                s[i] == '+' || s[i] == '-' || s[i] == '.' || s[i] == 'e' || s[i] == 'E')) ++i;
        if (i == j) { ok = false; return 0; }
        return std::stod(s.substr(j, i - j));
    }
    Tri triple() {
        if (!eat('[')) { ok = false; return {}; }
        int e = (int)num(); eat(',');
        double a = num(); eat(',');
        double b = num();
        if (!eat(']')) ok = false;
        return { e, a, b };
    }
    Tris arr() {
        Tris v;
        if (!eat('[')) { ok = false; return v; }
        if (eat(']')) return v;
        do { v.push_back(triple()); } while (eat(','));
        if (!eat(']')) ok = false;
        return v;
    }
    std::map<int, Tris> top() {
        std::map<int, Tris> m;
        if (!eat('{')) { ok = false; return m; }
        if (eat('}')) return m;
        do { std::string k = str(); if (!eat(':')) { ok = false; break; } m[std::stoi(k)] = arr(); } while (eat(','));
        if (!eat('}')) ok = false;
        return m;
    }
};

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: test_allowed_range <gridscale> <vectorfield> <mesh.obj> <fixture.json>\n"; return 2; }
    const double    scale = std::stod(argv[1]);
    const FieldType ft    = parse_field_type(argv[2]);
    const char*     mesh  = argv[3];
    const char*     fxp   = argv[4];

    std::ifstream f(fxp);
    CHECK(f.good());                                 // fixture opens
    std::stringstream ss; ss << f.rdbuf();
    std::string text = ss.str();
    JParser jp{ text };
    std::map<int, Tris> expected = jp.top();
    CHECK(jp.ok);                                    // fixture parses

    TmeshPipeline P(mesh, scale, ft);
    CHECK(P.ok);                                     // mesh loads + pipeline builds
    auto& tmm = *P.tmm;

    const double EPS = 1e-9;
    int n_ok = 0, n_fail = 0;
    for (auto& [tqid, exp_in] : expected) {
        CHECK(tqid >= 0 && tqid < (int)tmm.tquads.size());

        Tris got = tmm.allowed_range_tquads({tqid});
        Tris exp = exp_in;
        std::sort(got.begin(), got.end());
        std::sort(exp.begin(), exp.end());

        bool same = got.size() == exp.size();
        for (size_t k = 0; same && k < got.size(); ++k) {
            auto& [ge, gr0, gr1] = got[k];
            auto& [ee, er0, er1] = exp[k];
            if (ge != ee || std::abs(gr0 - er0) > EPS || std::abs(gr1 - er1) > EPS) same = false;
        }

        if (!same) {
            std::cerr << "FAIL: allowed_range regression  " << mesh << " tq" << tqid
                      << "  expected=" << exp.size() << " got=" << got.size() << "\n";
            // list eids present on only one side (quick diff)
            std::set<int> ge, ee;
            for (auto& [e, a, b] : got) ge.insert(e);
            for (auto& [e, a, b] : exp) ee.insert(e);
            auto report = [](const char* tag, const std::set<int>& a, const std::set<int>& b) {
                std::cerr << "  " << tag << ":"; int c = 0;
                for (int e : a) if (!b.count(e)) { std::cerr << " " << e; if (++c >= 40) { std::cerr << " ..."; break; } }
                std::cerr << "\n";
            };
            report("missing (in fixture, not produced)", ee, ge);
            report("extra   (produced, not in fixture)", ge, ee);
            ++n_fail;
        } else {
            ++n_ok;
        }
    }

    std::cout << "[test_allowed_range] " << mesh << "  ok=" << n_ok << " fail=" << n_fail << "\n";
    return n_fail ? 1 : 0;
}
