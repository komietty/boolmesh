//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::cast_abs_to_unsigned)]
#![allow(unused_braces)]

mod manifold;
mod triangulation;
mod simplification;
mod common;
mod boolean03;
mod boolean45;
mod compose;
mod tests;
mod cleanup;
use crate::boolean03::boolean03;
use crate::boolean45::boolean45;
use crate::simplification::simplify_topology;
use crate::triangulation::triangulate;
use crate::common::*;
use crate::cleanup::*;
use crate::manifold::*;

pub use crate::common::{Real, Vec2, Vec3, Vec4, Mat3, K_PRECISION};

pub trait Data: Clone + Send + Sync + std::fmt::Debug + PartialEq {}
impl<T> Data for T where T: Clone + Send + Sync + std::fmt::Debug + PartialEq {}

pub mod prelude {
    pub use crate::common::OpType;
    pub use crate::manifold::Manifold;
    pub use crate::compute_boolean;
    pub use crate::compose::{
        compose,
        fractal,
        extrude,
        generate_cone,
        generate_cube,
        generate_torus,
        generate_cylinder,
        generate_uv_sphere,
        generate_icosphere,
    };
}

pub fn compute_boolean<T: Data>(
    mp: &Manifold<T>,
    mq: &Manifold<T>,
    op: OpType,
) -> Result<Manifold<T>, String> {
    let eps = mp.eps.max(mq.eps);
    let tol = mp.tol.max(mq.tol);

    let     b03 = boolean03(mp, mq, &op);
    let mut b45 = boolean45(mp, mq, &b03, &op);
    let mut trg = triangulate(mp, mq, &b45, eps)?;

    simplify_topology(
        &mut trg.hs,
        &mut b45.ps,
        &mut trg.ns,
        &mut trg.rs,
        b45.nv_from_p,
        b45.nv_from_q,
        eps
    );

    cleanup_unused_verts(
        &mut b45.ps,
        &mut trg.hs,
        &mut trg.rs
    );

    let mut var = None;
    if let (Some(vp), Some(vq)) = (&mp.variable, &mq.variable) {
        var = Some(
            trg.rs.iter().map(|r| {
                if r.mid == 0 { vp[r.fid].clone() }
                else          { vq[r.fid].clone() }
            }).collect());
    }

    Manifold::new(
        b45.ps,
        trg.hs.iter().map(|h| h.tail).collect::<Vec<_>>(),
        var,
        Some(eps),
        Some(tol)
    )
}




