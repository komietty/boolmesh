//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

pub mod hmesh;
pub mod bounds;
pub mod collider;

use std::cmp::Ordering;
use std::fmt::Debug;
use bounds::BBox;
use hmesh::Hmesh;
use collider::{morton_code, MortonCollider};
use crate::{Data, Mat3, Half, Real, Vec3, Vec3u, K_PRECISION, next_of};
#[cfg(feature = "rayon")] use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct Manifold<T> {
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
    pub variable: Option<Vec<T>>, // variable on the mfd
    pub collider: MortonCollider, //
    pub coplanar: Vec<i32>,       // indices of coplanar faces
}

impl<T: Data> Manifold<T> {
    pub fn new<P: IntoVec3f, I: IntoVec3u>(
        pos: P,
        idx: I,
        var: Option<Vec<T>>,
        eps: Option<Real>,
        tol: Option<Real>,
    ) -> Result<Self, String> {
        let ps = pos.into_vec3f();
        let ts = idx.into_vec3u();
        let bb = BBox::new(None, &ps);
        let mut var = var;
        let mut mrt = compute_face_morton(&ps, &ts, &bb);
        let hm = sort_faces(&ps, &ts, &mut var, &mut mrt.0, &mut mrt.1)?;
        let hs = hm.half.iter().map(|&i| Half::new(hm.tail[i], hm.head[i], hm.twin[i])).collect::<Vec<_>>();

        let mut e = K_PRECISION * bb.scale();
        e = if e.is_finite() { e } else { -1. };
        let eps = if let Some(e_) = eps { e_ } else { e };
        let tol = if let Some(t_) = tol { t_ } else { e };
        let collider = MortonCollider::new(&mrt.0, &mrt.1);
        let coplanar = compute_coplanar_idx(&ps, &hm.fns, &hs, eps);

        let mfd = Manifold {
            ps,
            hs,
            eps,
            tol,
            nv: hm.nv,
            nf: hm.nf,
            nh: hm.nh,
            bounding_box: bb,
            vert_normals: hm.vns,
            face_normals: hm.fns,
            variable: var,
            collider,
            coplanar,
        };

        if !is_manifold(&mfd) { return Err("The input mesh is not manifold".into()); }
        Ok(mfd)
    }

    //fn set_epsilon(&mut self, min_epsilon: Real, use_single: bool) {
    //    let scl = self.bounding_box.scale();
    //    let mut e = min_epsilon.max(K_PRECISION * scl);
    //    e = if e.is_finite() { e } else { -1. };
    //    let t = if use_single { e.max(Real::EPSILON * scl) } else { e };
    //    self.eps = e;
    //    self.tol = self.tol.max(t);
    //}

    pub fn translate(&mut self, x: f64, y: f64, z: f64) {
        let t = Vec3::new(x as Real, y as Real, z as Real);
        *self = Manifold::new(
            self.ps.iter().map(|p| *p + t).collect::<Vec<_>>(),
            self.hs.iter().map(|h| h.tail).collect::<Vec<_>>(),
            self.variable.clone(),
            None,
            None
        ).unwrap();
    }

    pub fn rotate(&mut self, x: f64, y: f64, z: f64) {
        let r = Mat3::from_euler(glam::EulerRot::XYZ, x as Real, y as Real, z as Real);
        *self = Manifold::new(
            self.ps.iter().map(|p| r * *p).collect::<Vec<_>>(),
            self.hs.iter().map(|h| h.tail).collect::<Vec<_>>(),
            self.variable.clone(),
            None,
            None
        ).unwrap();
    }

    pub fn scale(&mut self, x: f64, y: f64, z: f64) {
        let s = Vec3::new(x as Real, y as Real, z as Real);
        *self = Manifold::new(
            self.ps.iter().map(|p| p * s).collect::<Vec<_>>(),
            self.hs.iter().map(|h| h.tail).collect::<Vec<_>>(),
            self.variable.clone(),
            None,
            None
        ).unwrap();
    }

    pub fn set_variable(&mut self, var: Vec<T>) {
        self.variable = Some(var);
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

fn sort_faces<T: Data>(
    pos: &[Vec3],
    idx: &[Vec3u],
    var: &mut Option<Vec<T>>,
    face_bboxes: &mut Vec<BBox>,
    face_morton: &mut Vec<u32>
) -> Result<Hmesh, String> {
    let mut map = (0..face_morton.len()).collect::<Vec<_>>();
    map.sort_by_key(|&i| face_morton[i]);

    *face_bboxes = map.iter().map(|&i| face_bboxes[i].clone()).collect::<Vec<_>>();
    *face_morton = map.iter().map(|&i| face_morton[i]).collect::<Vec<_>>();

    if let Some(v) = var {
        *v = map.iter().map(|&i| v[i].clone()).collect();
    }

    Hmesh::new(pos, &map.iter().map(|&i| idx[i]).collect::<Vec<_>>())
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

fn is_manifold<T: Data>(mfd: &Manifold<T>) -> bool {
    mfd.hs.iter().enumerate().all(|(i, h)| {
        if h.tail().is_none() || h.head().is_none() { return true; }
        match h.pair() {
            None => { false },
            Some(pair) => {
                let mut good = true;
                good &= mfd.hs[pair].pair() == Some(i);
                good &= h.tail != h.head;
                good &= h.tail == mfd.hs[pair].head;
                good &= h.head == mfd.hs[pair].tail;
                good
            }
        }
    })
}



pub trait IntoVec3f { fn into_vec3f(self) -> Vec<Vec3>;  }
pub trait IntoVec3u { fn into_vec3u(self) -> Vec<Vec3u>; }

impl IntoVec3f for &[f32]          { fn into_vec3f(self) -> Vec<Vec3> { self.chunks(3).map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &[f64]          { fn into_vec3f(self) -> Vec<Vec3> { self.chunks(3).map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &[[f32; 3]]     { fn into_vec3f(self) -> Vec<Vec3> { self.iter().map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &[[f64; 3]]     { fn into_vec3f(self) -> Vec<Vec3> { self.iter().map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &Vec<[f32; 3]>  { fn into_vec3f(self) -> Vec<Vec3> { self.iter().map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &Vec<[f64; 3]>  { fn into_vec3f(self) -> Vec<Vec3> { self.iter().map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &Vec<f32>       { fn into_vec3f(self) -> Vec<Vec3> { self.chunks(3).map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &Vec<f64>       { fn into_vec3f(self) -> Vec<Vec3> { self.chunks(3).map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for &[Vec3]         { fn into_vec3f(self) -> Vec<Vec3> { self.to_vec() } }
impl IntoVec3f for Vec<[f32; 3]>   { fn into_vec3f(self) -> Vec<Vec3> { self.iter().map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for Vec<[f64; 3]>   { fn into_vec3f(self) -> Vec<Vec3> { self.iter().map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for Vec<f32>        { fn into_vec3f(self) -> Vec<Vec3> { self.chunks(3).map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for Vec<f64>        { fn into_vec3f(self) -> Vec<Vec3> { self.chunks(3).map(|p| Vec3::new(p[0] as Real, p[1] as Real, p[2] as Real)).collect() }}
impl IntoVec3f for Vec<Vec3>       { fn into_vec3f(self) -> Vec<Vec3> { self }}

impl IntoVec3u for &[u32]           { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for &[u64]           { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for &[usize]         { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0], p[1], p[2])).collect() }}
impl IntoVec3u for &[[u32; 3]]      { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for &[[u64; 3]]      { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for &[[usize; 3]]    { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0], p[1], p[2])).collect() } }
impl IntoVec3u for &Vec<[u32; 3]>   { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for &Vec<[u64; 3]>   { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for &Vec<[usize; 3]> { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0], p[1], p[2])).collect() } }
impl IntoVec3u for &Vec<u32>        { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() } }
impl IntoVec3u for &Vec<u64>        { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() } }
impl IntoVec3u for &Vec<usize>      { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0], p[1], p[2])).collect() } }
impl IntoVec3u for &[Vec3u]         { fn into_vec3u(self) -> Vec<Vec3u> { self.to_vec() } }
impl IntoVec3u for Vec<[u32; 3]>    { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for Vec<[u64; 3]>    { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() }}
impl IntoVec3u for Vec<[usize; 3]>  { fn into_vec3u(self) -> Vec<Vec3u> { self.iter().map(|p| Vec3u::new(p[0], p[1], p[2])).collect() } }
impl IntoVec3u for Vec<u32>         { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() } }
impl IntoVec3u for Vec<u64>         { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0] as usize, p[1] as usize, p[2] as usize)).collect() } }
impl IntoVec3u for Vec<usize>       { fn into_vec3u(self) -> Vec<Vec3u> { self.chunks(3).map(|p| Vec3u::new(p[0], p[1], p[2])).collect() } }
impl IntoVec3u for Vec<Vec3u>       { fn into_vec3u(self) -> Vec<Vec3u> { self } }

