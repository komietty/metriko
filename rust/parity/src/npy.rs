//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Reader of the .npy files written by tools/dump: little endian `<f8`, `<i4`, `<c16`, `|u1`, C order.

use metriko::Complex;
use std::path::Path;

#[derive(Clone, Debug)]
pub enum NpyData {
    F8(Vec<f64>),
    I4(Vec<i32>),
    C16(Vec<Complex>),
    U1(Vec<u8>),
}

#[derive(Clone, Debug)]
pub struct Npy {
    pub name: String,
    pub shape: Vec<usize>,
    pub data: NpyData,
}

impl Npy {
    pub fn read(path: &Path) -> Result<Npy, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::parse(&bytes).map(|(shape, data)| Npy { name: path.file_name().unwrap().to_string_lossy().into_owned(), shape, data })
            .map_err(|e| format!("{}: {e}", path.display()))
    }

    fn parse(b: &[u8]) -> Result<(Vec<usize>, NpyData), String> {
        if b.len() < 10 || &b[..6] != b"\x93NUMPY" { return Err("not an npy file".into()); }
        let (hl, off) = match b[6] {
            1 => (u16::from_le_bytes([b[8], b[9]]) as usize, 10),
            2 | 3 => (u32::from_le_bytes([b[8], b[9], b[10], b[11]]) as usize, 12),
            v => return Err(format!("npy version {v} not supported")),
        };
        let header = std::str::from_utf8(&b[off..off + hl]).map_err(|e| e.to_string())?;
        let field = |key: &str| -> Result<&str, String> {
            let i = header.find(&format!("'{key}'")).ok_or(format!("no {key} in header {header}"))?;
            Ok(header[i + key.len() + 2..].trim_start().trim_start_matches(':').trim_start())
        };
        let descr = field("descr")?;
        let descr = &descr[1..descr[1..].find('\'').ok_or("bad descr")? + 1];
        if field("fortran_order")?.starts_with("True") { return Err("fortran order not supported".into()); }
        let shape = field("shape")?;
        let shape = &shape[1..shape.find(')').ok_or("bad shape")?];
        let shape: Vec<usize> = shape.split(',').map(str::trim).filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>().map_err(|e| e.to_string())).collect::<Result<_, _>>()?;
        let n: usize = shape.iter().product();
        let raw = &b[off + hl..];
        let need = |size: usize| if raw.len() < n * size { Err(format!("{} bytes of data, {} expected", raw.len(), n * size)) } else { Ok(()) };
        let data = match descr {
            "<f8" => { need(8)?; NpyData::F8(raw.chunks_exact(8).take(n).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect()) }
            "<i4" => { need(4)?; NpyData::I4(raw.chunks_exact(4).take(n).map(|c| i32::from_le_bytes(c.try_into().unwrap())).collect()) }
            "<c16" => {
                need(16)?;
                NpyData::C16(raw.chunks_exact(16).take(n).map(|c| Complex::new(
                    f64::from_le_bytes(c[..8].try_into().unwrap()),
                    f64::from_le_bytes(c[8..].try_into().unwrap()))).collect())
            }
            "|u1" | "<u1" => { need(1)?; NpyData::U1(raw[..n].to_vec()) }
            d => return Err(format!("dtype {d} not supported")),
        };
        Ok((shape, data))
    }

    pub fn len(&self) -> usize { self.shape.iter().product() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }
    pub fn rows(&self) -> usize { self.shape.first().copied().unwrap_or(1) }
    /// columns of a 2D array, 1 for a vector
    pub fn cols(&self) -> usize { if self.shape.len() >= 2 { self.shape[1..].iter().product() } else { 1 } }

    pub fn f64s(&self) -> &[f64] {
        match &self.data { NpyData::F8(v) => v, d => panic!("{}: f8 expected, found {}", self.name, kind(d)) }
    }
    pub fn i32s(&self) -> &[i32] {
        match &self.data { NpyData::I4(v) => v, d => panic!("{}: i4 expected, found {}", self.name, kind(d)) }
    }
    pub fn c64s(&self) -> &[Complex] {
        match &self.data { NpyData::C16(v) => v, d => panic!("{}: c16 expected, found {}", self.name, kind(d)) }
    }
    pub fn u8s(&self) -> &[u8] {
        match &self.data { NpyData::U1(v) => v, d => panic!("{}: u1 expected, found {}", self.name, kind(d)) }
    }
    pub fn bools(&self) -> Vec<bool> { self.u8s().iter().map(|&b| b != 0).collect() }

    /// a 2D f8 array as rows of N
    pub fn rows_f64<const N: usize>(&self) -> Vec<[f64; N]> {
        assert_eq!(self.cols(), N, "{}: {} columns, {N} expected", self.name, self.cols());
        self.f64s().chunks_exact(N).map(|c| c.try_into().unwrap()).collect()
    }
    /// a 2D i4 array as rows of N
    pub fn rows_i32<const N: usize>(&self) -> Vec<[i32; N]> {
        assert_eq!(self.cols(), N, "{}: {} columns, {N} expected", self.name, self.cols());
        self.i32s().chunks_exact(N).map(|c| c.try_into().unwrap()).collect()
    }
    /// an (n, 3) f8 array as points
    pub fn dvec3s(&self) -> Vec<glam::DVec3> {
        self.rows_f64::<3>().into_iter().map(glam::DVec3::from_array).collect()
    }
}

fn kind(d: &NpyData) -> &'static str {
    match d { NpyData::F8(_) => "f8", NpyData::I4(_) => "i4", NpyData::C16(_) => "c16", NpyData::U1(_) => "u1" }
}
