//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
#![allow(clippy::unnecessary_cast)]

pub mod cone;
pub mod cube;
pub mod cylinder;
pub mod extrusion;
pub mod sphere;
pub mod torus;

pub use cone::*;
pub use cube::*;
pub use cylinder::*;
pub use extrusion::*;
pub use sphere::*;
pub use torus::*;

use crate::manifold::ManifoldError;
use crate::Manifold;

pub fn compose(ms: &Vec<Manifold>) -> Result<Manifold, ManifoldError> {
    let mut ps = vec![];
    let mut ts = vec![];
    let mut offset = 0;
    for m in ms {
        for h in m.hs.iter() {
            ts.push(h.tail + offset);
        }
        for p in m.ps.iter() {
            ps.push(p.x as f64);
            ps.push(p.y as f64);
            ps.push(p.z as f64);
        }
        offset += m.nv;
    }
    Manifold::new(&ps, &ts)
}

pub fn fractal(
    hole: &Manifold,
    holes: &mut Vec<Manifold>,
    x: f64,
    y: f64,
    w: f64,
    depth: usize,
    depth_max: usize,
) {
    let w = w / 3.;
    let mut m = hole.clone();
    m.scale(w, w, 1.);
    m.translate(x, y, 0.);
    holes.push(m);

    if depth == depth_max {
        return;
    }

    for xy in [
        (x - w, y - w),
        (x - w, y),
        (x - w, y + w),
        (x, y + w),
        (x + w, y + w),
        (x + w, y),
        (x + w, y - w),
        (x, y - w),
    ] {
        fractal(hole, holes, xy.0, xy.1, w, depth + 1, depth_max);
    }
}
