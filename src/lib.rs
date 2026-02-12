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

use crate::boolean03::boolean03;
use crate::boolean45::boolean45;
use crate::simplification::simplify_topology;
use crate::triangulation::triangulate;
use crate::common::*;
use crate::manifold::*;

pub use crate::common::{Real, Vec2, Vec3, Vec4, Mat3, K_PRECISION};
use crate::manifold::bounds::BBox;
use crate::manifold::collider::{morton_code, K_NO_CODE};

pub mod prelude {
    pub use crate::common::OpType;
    pub use crate::manifold::Manifold;
    pub use crate::compute_boolean;
    pub use crate::compute_boolean_with_attributes;
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

//==================== simple boolean ====================//

pub fn compute_boolean(
    mp: &Manifold,
    mq: &Manifold,
    op: OpType,
) -> Result<Manifold, String> {
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

    Manifold::new_impl(
        b45.ps,
        trg.hs.chunks(3).map(|hs| Vec3u::new(hs[0].tail, hs[1].tail, hs[2].tail)).collect(),
        Some(eps),
        Some(tol)
    )
}

//==================== attribute boolean ====================//

pub trait Attr: Clone + Send + Sync + std::fmt::Debug + PartialEq {}
impl<T> Attr for T where T: Clone + Send + Sync + std::fmt::Debug + PartialEq {}

#[derive(Clone, Debug)]
pub struct AttrManifold<T> {
    pub manifold: Manifold,
    pub attribute: Vec<T>,
}

impl<T: Attr> AttrManifold<T> {
    pub fn new(mfd: Manifold, attr: Vec<T>, sort: bool) -> Self {
        if sort {
            let attr = mfd.sort_map.iter().map(|&i| attr[i].clone()).collect::<Vec<_>>();
            return Self { manifold: mfd, attribute: attr };
        }
        Self { manifold: mfd, attribute: attr }
    }
}

pub fn compute_boolean_with_attributes<T: Attr>(
    attr_mp: &AttrManifold<T>,
    attr_mq: &AttrManifold<T>,
    op     : OpType,
) -> Result<AttrManifold<T>, String> {
    let eps = attr_mp.manifold.eps.max(attr_mq.manifold.eps);
    let tol = attr_mp.manifold.tol.max(attr_mq.manifold.tol);
    let mp = &attr_mp.manifold;
    let mq = &attr_mq.manifold;

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

    let mfd = Manifold::new_impl(
        b45.ps,
        trg.hs.chunks(3).map(|hs| Vec3u::new(hs[0].tail, hs[1].tail, hs[2].tail)).collect(),
        Some(eps),
        Some(tol)
    )?;

    let mut attr = vec![];
    let ap = &attr_mp.attribute;
    let aq = &attr_mq.attribute;
    if !ap.is_empty() && !aq.is_empty() {
        attr = trg.rs.iter().map(|r| {
            if r.mid == 0 { ap[r.fid].clone() }
            else          { aq[r.fid].clone() }
        }).collect();
    }

    Ok(AttrManifold::new(mfd, attr, true))
}

//==================== functions ====================//

fn cleanup_unused_verts(
    ps: &mut Vec<Vec3>,
    hs: &mut Vec<Half>,
    rs: &mut Vec<Tref>,
) {
    let bb = BBox::new(None, ps);
    let mt = ps.iter().map(|p| morton_code(p, &bb)).collect::<Vec<_>>();

    let mut new2old = (0..ps.len()).collect::<Vec<_>>();
    let mut old2new = vec![0; ps.len()];
    new2old.sort_by_key(|&i| mt[i]);
    for (new, &old) in new2old.iter().enumerate() { old2new[old] = new; }

    // reindex verts
    for h in hs.iter_mut() {
        if h.pair().is_none() { continue; }
        h.tail = old2new[h.tail];
        h.head = old2new[h.head];
    }

    // truncate pos container
    let nv = new2old
        .iter()
        .position(|&v| mt[v] >= K_NO_CODE)
        .unwrap_or(new2old.len());

    new2old.truncate(nv);

    *ps = new2old.iter().map(|&i| ps[i]).collect();
    *rs = hs.chunks(3).enumerate().filter_map(|(i, t)| t[0].pair().map(|_| rs[i])).collect();
    *hs = hs.iter().filter(|h| h.pair().is_some()).cloned().collect();
}
