//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stage 01: `Hmesh(V / unit, F)` with its geometry, bit for bit.

use parity::*;
use parity::stage::*;

#[test]
fn hmesh_with_geometry() {
    for (c, g) in golden_cases() {
        let hm = hmesh(&g, "01_hm", true);
        check_hmesh(&c.name(), &hm, &g, "01_hm", true).unwrap_or_else(|e| panic!("{e}"));
    }
}
