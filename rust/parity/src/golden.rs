//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Discovery of the golden dumps: testdata/cases.json lists the cases, tools/golden.sh writes them to
//! testdata/golden/<first 7 chars of rust/PORTED_FROM>/<stem>@<scale>/.

use crate::npy::Npy;
use serde_json::Value;
use std::path::{Path, PathBuf};

/// the metriko repository (the C++ root, parent of rust/); $METRIKO_ROOT overrides it
pub fn repo_root() -> PathBuf {
    match std::env::var_os("METRIKO_ROOT") {
        Some(d) => PathBuf::from(d),
        None => Path::new(env!("CARGO_MANIFEST_DIR")).join("../..").canonicalize().expect("repository root"),
    }
}

/// the C++ commit the Rust tree mirrors (rust/PORTED_FROM)
pub fn ported_from() -> String {
    let p = repo_root().join("rust/PORTED_FROM");
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).trim().to_string()
}

fn cases_json() -> Value {
    let p = repo_root().join("testdata/cases.json");
    serde_json::from_str(&std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display())))
        .unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

/// $METRIKO_MODELS, or the default of testdata/cases.json
pub fn models_dir() -> PathBuf {
    match std::env::var_os("METRIKO_MODELS") {
        Some(d) => PathBuf::from(d),
        None => PathBuf::from(cases_json()["models_default"].as_str().expect("models_default")),
    }
}

/// $METRIKO_GOLDEN, or testdata/golden/<short commit>
pub fn golden_root() -> PathBuf {
    match std::env::var_os("METRIKO_GOLDEN") {
        Some(d) => PathBuf::from(d),
        None => repo_root().join("testdata/golden").join(&ported_from()[..7]),
    }
}

#[derive(Clone, Debug)]
pub struct Case {
    pub mesh: String,      // file name in the models directory
    pub stem: String,
    pub scale_str: String, // as written in cases.json, also the directory suffix
    pub scale: f64,        // parsed as the C++ std::stod does (both are correctly rounded)
}

impl Case {
    pub fn name(&self) -> String { format!("{}@{}", self.stem, self.scale_str) }
    pub fn mesh_path(&self) -> PathBuf { models_dir().join(&self.mesh) }
    pub fn golden_dir(&self) -> PathBuf { golden_root().join(self.name()) }
    pub fn golden(&self) -> Option<Golden> { Golden::open(&self.golden_dir()) }
}

/// every case of testdata/cases.json
pub fn cases() -> Vec<Case> {
    cases_json()["cases"].as_array().expect("cases").iter().map(|c| {
        let mesh = c["mesh"].as_str().expect("mesh").to_string();
        let scale_str = c["scale"].as_str().expect("scale must be a string").to_string();
        Case {
            stem: Path::new(&mesh).file_stem().unwrap().to_string_lossy().into_owned(),
            scale: scale_str.parse().expect("scale"),
            mesh,
            scale_str,
        }
    }).collect()
}

/// the cases whose golden dump exists (and whose dump completed: summary.json present). Prints what is missing, and
/// panics when there is none at all, so that a test never passes vacuously; set METRIKO_ALLOW_NO_GOLDEN=1 to allow it
pub fn golden_cases() -> Vec<(Case, Golden)> {
    let mut res = vec![];
    for c in cases() {
        match c.golden() {
            Some(g) => res.push((c, g)),
            None => eprintln!("[parity] no golden dump for {} in {} (run tools/golden.sh)", c.name(), golden_root().display()),
        }
    }
    if res.is_empty() && std::env::var_os("METRIKO_ALLOW_NO_GOLDEN").is_none() {
        panic!("no golden dump found under {}: run tools/golden.sh", golden_root().display());
    }
    res
}

/// one case directory written by tools/dump
pub struct Golden {
    pub dir: PathBuf,
    pub meta: serde_json::Map<String, Value>,
    pub summary: Value,
}

impl Golden {
    pub fn open(dir: &Path) -> Option<Golden> {
        let summary = std::fs::read_to_string(dir.join("summary.json")).ok()?;
        let manifest: Value = serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).ok()?).ok()?;
        Some(Golden {
            dir: dir.to_path_buf(),
            meta: manifest["meta"].as_object().cloned().unwrap_or_default(),
            summary: serde_json::from_str(&summary).ok()?,
        })
    }

    /// `key` is "<stage>.<name>", e.g. "02_rawf.field" for 02_rawf.field.npy
    pub fn load(&self, key: &str) -> Npy {
        Npy::read(&self.dir.join(format!("{key}.npy"))).unwrap_or_else(|e| panic!("{e}"))
    }

    pub fn has(&self, key: &str) -> bool { self.dir.join(format!("{key}.npy")).exists() }

    /// a manifest meta entry, e.g. "00_input.unit"
    pub fn meta(&self, key: &str) -> &Value {
        self.meta.get(key).unwrap_or_else(|| panic!("{}: no meta {key}", self.dir.display()))
    }

    pub fn meta_f64(&self, key: &str) -> f64 { self.meta(key).as_f64().unwrap_or_else(|| panic!("meta {key} is not a number")) }
    pub fn meta_i64(&self, key: &str) -> i64 { self.meta(key).as_i64().unwrap_or_else(|| panic!("meta {key} is not an integer")) }
}
