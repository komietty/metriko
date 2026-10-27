//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of metriko_web/wasm/cli.cpp.

//! metriko <mesh.obj> <scale> [positions.bin]
//!
//! prints the same summary as metriko_web's native cli (and wasm/test.mjs), so that both can be compared line by
//! line; writes the quad vertices as raw little endian float64 xyz if a file is given.

use std::io::Write;
use std::path::Path;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: metriko <mesh> <scale> [positions.bin]");
        std::process::exit(2);
    }
    let m = match parity::load_mesh(Path::new(&args[1])) {
        Ok(m) => m,
        Err(e) => { eprintln!("failed: {e}"); std::process::exit(1); }
    };
    let scale: f64 = match args[2].parse() {
        Ok(s) => s,
        Err(e) => { eprintln!("failed: scale {}: {e}", args[2]); std::process::exit(2); }
    };

    let t0 = Instant::now();
    let res = metriko::compute_quadrangulation(&m.v, &m.f, scale);
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    let q = match res {
        Ok(q) => q,
        Err(e) => { eprintln!("failed: {e}"); std::process::exit(1); }
    };

    let s = parity::Summary::new(m.v.len(), m.f.len(), &q.pos, &q.idx, &q.val);
    for l in s.lines() { println!("{l}"); }
    println!("time    {ms:.0} ms");
    if let Some(out) = args.get(3) {
        let bytes: Vec<u8> = q.pos.iter().flat_map(|p| p.to_array()).flat_map(f64::to_le_bytes).collect();
        if let Err(e) = std::fs::File::create(out).and_then(|mut f| f.write_all(&bytes)) {
            eprintln!("failed: {out}: {e}");
            std::process::exit(1);
        }
    }
}
