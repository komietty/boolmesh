//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::cast_abs_to_unsigned)]
#![allow(unused_braces)]

mod boolean03;
mod boolean45;
mod common;
mod compose;
mod manifold;
mod simplification;
mod tests;
mod triangulation;

use thiserror::Error;

use crate::boolean03::boolean03;
use crate::boolean45::boolean45;
use crate::common::*;
use crate::manifold::*;
use crate::simplification::simplify_topology;
use crate::triangulation::triangulate;
use crate::triangulation::TriangulationError;

pub use crate::common::{Mat3, Real, Vec2, Vec3, Vec4, K_PRECISION};

pub mod prelude {
    pub use crate::common::OpType;
    pub use crate::compose::{
        compose, extrude, fractal, generate_cone, generate_cube, generate_cylinder,
        generate_icosphere, generate_torus, generate_uv_sphere,
    };
    pub use crate::compute_boolean;
    pub use crate::manifold::Manifold;
}

pub fn compute_boolean(mp: &Manifold, mq: &Manifold, op: OpType) -> Result<Manifold, BooleanError> {
    let eps = mp.eps.max(mq.eps);
    let tol = mp.tol.max(mq.tol);

    let b03 = boolean03(mp, mq, &op);
    let mut b45 = boolean45(mp, mq, &b03, &op);
    let mut trg = triangulate(mp, mq, &b45, eps)?;

    simplify_topology(
        &mut trg.hs,
        &mut b45.ps,
        &mut trg.ns,
        &mut trg.rs,
        b45.nv_from_p,
        b45.nv_from_q,
        eps,
    );

    cleanup_unused_verts(&mut b45.ps, &mut trg.hs);

    let manifold = Manifold::new_impl(
        b45.ps,
        trg.hs
            .chunks(3)
            .map(|hs| Vec3u::new(hs[0].tail, hs[1].tail, hs[2].tail))
            .collect(),
        Some(eps),
        Some(tol),
    )?;

    Ok(manifold)
}

#[derive(Debug, Error)]
pub enum BooleanError {
    #[error("{0}")]
    Trangulate(#[from] TriangulationError),

    #[error("{0}")]
    Manifold(#[from] ManifoldError),
}

//pub fn compute_boolean_from_raw_data(
//    pos0: &[Real],
//    idx0: &[usize],
//    pos1: &[Real],
//    idx1: &[usize],
//    op_type: usize
//) -> Result<Manifold, String>{
//    let mp = Manifold::new(&pos0, &idx0)?;
//    let mq = Manifold::new(&pos1, &idx1)?;
//    let op = match op_type {
//        0 => OpType::Add,
//        1 => OpType::Subtract,
//        2 => OpType::Intersect,
//        _ => return Err("Invalid op_type".into())
//    };
//    compute_boolean(&mp, &mq, op)
//}
