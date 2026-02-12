//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

pub mod hmesh;
pub mod bounds;
pub mod collider;

use std::cmp::Ordering;
use std::collections::HashMap;
use bounds::BBox;
use crate::collider::{morton_code, MortonCollider, K_NO_CODE};
use crate::{Real, Half, Vec3, Vec3u, K_PRECISION, next_of, Mat3};
use super::hmesh::Hmesh;
#[cfg(feature = "rayon")] use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct Manifold {
    pub ps: Vec<Vec3>,            // positions
    pub hs: Vec<Half>,            // halfedges
    pub nv: usize,                // number of vertices
    pub nf: usize,                // number of faces
    pub nh: usize,                // number of halfedges
    pub eps: Real,                // epsilon
    pub tol: Real,                // tolerance
    pub bounding_box: BBox,       //
    pub face_normals: Vec<Vec3>,  //
    pub vert_normals: Vec<Vec3>,  //
    pub sort_map: Vec<usize>,     //
    pub collider: MortonCollider, //
    pub coplanar: Vec<i32>,       // indices of coplanar faces
}

impl Manifold {
    pub fn new(pos: &[f64], idx: &[usize]) -> Result<Self, String> {

        if pos.len() % 3 != 0 { return Err("pos must be a multiple of 3".into()); }
        if idx.len() % 3 != 0 { return Err("idx must be a multiple of 3".into()); }

        // dedup vertices
        let mut hash  = HashMap::with_capacity(pos.len() / 3);
        let mut weld  = Vec::with_capacity(pos.len() / 3);
        let mut rmap = vec![0; pos.len()];

        for (i, p) in pos.chunks(3).enumerate() {
            let v = Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real);
            let k = (v.x.to_bits(), v.y.to_bits(), v.z.to_bits());
            if let Some(&w) = hash.get(&k) { rmap[i] = w; }
            else {
                let n = weld.len();
                weld.push(v);
                hash.insert(k, n);
                rmap[i] = n;
            }
        }

        // remove collapsed triangles
        let idx = idx.chunks(3).filter_map(|i| {
            let v = Vec3u::new(rmap[i[0]], rmap[i[1]], rmap[i[2]]);
            (v.x != v.y && v.y != v.z && v.z != v.x).then_some(v)
        }).collect::<Vec<_>>();

        Self::new_impl(weld, idx, None, None)
    }

    pub fn new_impl(
        ps : Vec<Vec3>,
        idx: Vec<Vec3u>,
        eps: Option<Real>,
        tol: Option<Real>,
    ) -> Result<Self, String> {
        let bb = BBox::new(None, &ps);
        let (mut f_bb, mut f_mt) = compute_face_morton(&ps, &idx, &bb);

        // sort faces by morton code
        let mut map = (0..f_mt.len()).collect::<Vec<_>>();
        map.sort_by_key(|&i| f_mt[i]);
        f_bb = map.iter().map(|&i| f_bb[i].clone()).collect::<Vec<_>>();
        f_mt = map.iter().map(|&i| f_mt[i]).collect::<Vec<_>>();

        let hm = Hmesh::new(&ps, &map.iter().map(|&i| idx[i]).collect::<Vec<_>>())?;
        let hs = hm.half.iter().map(|&i| Half::new(hm.tail[i], hm.head[i], hm.twin[i])).collect::<Vec<_>>();

        let mut e = K_PRECISION * bb.scale();
        e = if e.is_finite() { e } else { -1. };
        let eps = if let Some(e_) = eps { e_ } else { e };
        let tol = if let Some(t_) = tol { t_ } else { e };
        let collider = MortonCollider::new(&f_bb, &f_mt);
        let coplanar = compute_coplanar_idx(&ps, &hm.fns, &hs, eps);

        let mfd = Manifold {
            nv: hm.nv,
            nf: hm.nf,
            nh: hm.nh,
            ps,
            hs,
            bounding_box: bb,
            vert_normals: hm.vns,
            face_normals: hm.fns,
            sort_map: map,
            collider,
            coplanar,
            eps,
            tol,
        };

        if !mfd.is_manifold() { return Err("The input mesh is not manifold".into()); }
        Ok(mfd)
    }

    pub fn get_indices(&self) -> Vec<Vec3u> {
        self.hs.chunks(3).map(|cs| Vec3u::new(cs[0].tail, cs[1].tail, cs[2].tail)).collect()
    }

    pub fn set_epsilon(&mut self, min_epsilon: Real, use_single: bool) {
        let scl = self.bounding_box.scale();
        let mut e = min_epsilon.max(K_PRECISION * scl);
        e = if e.is_finite() { e } else { -1. };
        let t = if use_single { e.max(Real::EPSILON * scl) } else { e };
        self.eps = e;
        self.tol = self.tol.max(t);
    }

    pub fn is_manifold(&self) -> bool {
        self.hs.iter().enumerate().all(|(i, h)| {
            if h.tail().is_none() || h.head().is_none() { return true; }
            match h.pair() {
                None => { false },
                Some(pair) => {
                    let mut good = true;
                    good &= self.hs[pair].pair() == Some(i);
                    good &= h.tail != h.head;
                    good &= h.tail == self.hs[pair].head;
                    good &= h.head == self.hs[pair].tail;
                    good
                }
            }
        })
    }

    pub fn translate(&mut self, x: f64, y: f64, z: f64) {
        let t = Vec3::new(x as Real, y as Real, z as Real);
        let p = self.ps.iter().map(|p| *p + t).collect();
        *self = Manifold::new_impl(p, self.get_indices(), None, None).unwrap();
    }

    pub fn rotate(&mut self, x: f64, y: f64, z: f64) {
        let r = Mat3::from_euler(glam::EulerRot::XYZ, x as Real, y as Real, z as Real);
        let p = self.ps.iter().map(|p| r * *p).collect();
        *self = Manifold::new_impl(p, self.get_indices(), None, None).unwrap();
    }

    pub fn scale(&mut self, x: f64, y: f64, z: f64) {
        let p = self.ps.iter().map(|p| Vec3::new(p.x * x as Real, p.y * y as Real, p.z * z as Real)).collect();
        *self = Manifold::new_impl(p, self.get_indices(), None, None).unwrap();
    }
}

fn compute_face_morton(
    pos: &[Vec3],
    idx: &[Vec3u],
    bb: &BBox
) -> (Vec<BBox>, Vec<u32>) {
    let n = idx.len();
    let mut bbs = vec![BBox::default(); n];
    let mut mts = vec![0; n];

    #[cfg(feature = "rayon")] {
        bbs.par_iter_mut()
            .zip(mts.par_iter_mut())
            .zip(idx.par_iter())
            .for_each(|((bb_, mt_), f)| {
                let p0 = pos[f.x];
                let p1 = pos[f.y];
                let p2 = pos[f.z];
                bb_.union(&p0);
                bb_.union(&p1);
                bb_.union(&p2);
                *mt_ = morton_code(&((p0 + p1 + p2) / 3.), bb);
            });
    }

    #[cfg(not(feature = "rayon"))] {
        for (i, f) in idx.iter().enumerate() {
            let p0 = pos[f.x];
            let p1 = pos[f.y];
            let p2 = pos[f.z];
            bbs[i].union(&p0);
            bbs[i].union(&p1);
            bbs[i].union(&p2);
            mts[i] = morton_code(&((p0 + p1 + p2) / 3.), bb);
        }
    }


    (bbs, mts)
}

fn compute_coplanar_idx(
    ps: &[Vec3],
    ns: &[Vec3],
    hs: &[Half],
    tol: Real
) -> Vec<i32> {
    let nt = hs.len() / 3;
    let mut priority = vec![];
    let mut res = vec![-1; nt];

    for t in 0..nt {
        let i = t * 3;
        let area = if hs[i].tail().is_none() { 0.} else {
            let p0 = ps[hs[i].tail];
            let p1 = ps[hs[i].head];
            let p2 = ps[hs[i + 1].head];
            (p1 - p0).cross(p2 - p0).length_squared()
        };
        priority.push((area, t));
    }

    priority.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));

    let mut interior = vec![];
    for (_, t) in priority.iter() {
        if res[*t] != -1 { continue; }
        res[*t] = *t as i32;

        let i = t * 3;
        let p = ps[hs[i].tail];
        let n = ns[*t];

        interior.clear();
        interior.extend_from_slice(&[i, i + 1, i + 2]);

        while let Some(hi) = interior.pop() {
            let h1 = next_of(hs[hi].pair);
            let t1 = h1 / 3;

            if res[t1] != -1 { continue; }

            if (ps[hs[h1].head] - p).dot(n).abs() < tol {
                res[t1] = *t as i32;
                if interior.last().copied() == Some(hs[h1].pair) { interior.pop(); }
                else { interior.push(h1); }
                interior.push(next_of(h1));
            }
        }
    }
    res
}
