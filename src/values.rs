//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use std::collections::HashMap;
use crate::{Data, Real, Vec3, Half, Tref, Manifold};

fn nxt(i: usize) -> usize { (i + 1) % 3 }
fn prv(i: usize) -> usize { (i + 2) % 3 }

fn get_barycentric_ratio(v: &Vec3, pos: &[Vec3; 3], tol: Real) -> Vec3 {
    let mut uvw = Vec3::ZERO;
    let dif = [pos[2] - pos[1], pos[0] - pos[2], pos[1] - pos[0]];
    let lsq = [dif[0].length_squared(), dif[1].length_squared(), dif[2].length_squared()];
    let crs = dif[0].cross(dif[1]);
    let long = if lsq[0] > lsq[1] && lsq[0] > lsq[2] { 0 } else if lsq[1] > lsq[2] { 1 } else { 2 };
    let tol2 = tol * tol;

    for i in 0..3 { if (*v - pos[i]).length_squared() < tol2 { uvw[i] = 1.; return uvw; } }

    if lsq[long] < tol2 { Vec3::new(1., 0., 0.) }
    else if crs.length_squared() > lsq[long] * tol2 {
        for i in 0..3 {
            let c = dif[i].cross(*v - pos[nxt(i)]);
            uvw[i] = if c.length_squared() < lsq[i] * tol2 { 0. } else { c.dot(crs) };
        }
        // Within `tol` of all three edge lines (a thin triangle at a coarse tolerance)
        // every weight is zero; dividing would make every value NaN. Take the nearest
        // corner instead, still a valid interpolation.
        let sum = uvw[0] + uvw[1] + uvw[2];
        if sum.abs() > Real::MIN_POSITIVE && sum.is_finite() {
            uvw /= sum;
            return uvw;
        }
        let d = [0, 1, 2].map(|i| (*v - pos[i]).length_squared());
        let near = if d[0] <= d[1] && d[0] <= d[2] { 0 } else if d[1] <= d[2] { 1 } else { 2 };
        uvw = Vec3::ZERO;
        uvw[near] = 1.;
        uvw
    } else {
        let next = nxt(long);
        let last = nxt(next);
        let r = (*v - pos[next]).dot(dif[long]) / lsq[long];
        uvw[long] = 0.;
        uvw[next] = 1. - r;
        uvw[last] = r;
        uvw
    }
}

pub fn annotate<S: Data>(
    hs: &mut [Half],
    ps: &[Vec3],
    rs: &[Tref],
    mp: &Manifold<S>,
    mq: &Manifold<S>,
    eps: Real,
) -> Vec<S> {
    if mp.vals.is_empty() && mq.vals.is_empty() { return vec![]; }

    let none = usize::MAX;
    let mut vs = vec![];
    let mut shared = HashMap::new();

    for tri in 0..hs.len() / 3 {
        if hs[3 * tri..3 * tri + 3].iter().any(|h| h.tail == none) { continue; }
        let r = rs[tri];
        let pq = r.mid == 0;
        let m = if pq { mp } else { mq };
        let old = |j: usize| m.hs[3 * r.fid + j].prop;
        let val = |j: usize| m.vals[3 * r.fid + j].clone();
        let pos = [0, 1, 2].map(|j| m.ps[m.hs[3 * r.fid + j].tail]);

        for i in 0..3 {
            let h = 3 * tri + i;
            let v = hs[h].tail;
            let r = get_barycentric_ratio(&ps[v], &pos, eps);
            let w = |j: usize| val(j) * r[j];
            let k = if m.vals.is_empty() { (pq, none, none, none) }
                    else if let Some(j) = (0..3).find(|&j| r[j] == 1.) { (pq, none, old(j), none) }
                    else if let Some(j) = (0..3).rev().find(|&j| r[j] == 0.) { let (a, b) = (old(nxt(j)), old(prv(j))); (pq, v, a.min(b), a.max(b)) }
                    else { (pq, v, none, none) };

            hs[h].prop = *shared.entry(k).or_insert_with(|| {
                vs.push(if m.vals.is_empty() { S::default() } else { w(0) + w(1) + w(2) });
                vs.len() - 1
            });
        }
    }
    vs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_near_all_three_edges_of_thin_triangle() {
        // height 1.05e-4 > tol, so the triangle is not treated as degenerate, yet the point
        // is within tol = 1e-4 of all three edge lines and of no corner
        let pos = [Vec3::new(0., 0., 0.), Vec3::new(1., 0., 0.), Vec3::new(0.5, 1.05e-4, 0.)];
        let uvw = get_barycentric_ratio(&Vec3::new(0.45, 0.5e-4, 0.), &pos, 1e-4);
        assert!(uvw.is_finite(), "{uvw:?}");
        assert!((uvw[0] + uvw[1] + uvw[2] - 1.).abs() < 1e-12, "{uvw:?}");
    }
}
