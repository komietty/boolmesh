//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#[cfg(test)]
mod test_intersection {
    use crate::boolean03::Boolean03;
    use crate::boolean03::kernel03::winding03;
    use crate::boolean03::kernel12::intersect12;
    use crate::boolean45::boolean45;
    use crate::{OpType, Vec3};
    use crate::Manifold;

    pub fn gen_tet_a() -> Manifold {
        Manifold::new(
            &vec![-0.866025, -1., 0.5, 0., -1., -1., 0.866025, -1., 0.5, 0., 1., 0.],
            &vec![0, 3, 1, 1, 2, 0, 1, 3, 2, 2, 3, 0],
        ).unwrap()
    }

    pub fn gen_tet_b() -> Manifold {
        Manifold::new(
            &vec![-1., -0.866025, 0.5, -1., 0., -1., -1., 0.866025, 0.5, 1., 0., 0.],
            &vec![1, 3, 0, 1, 0, 2, 2, 3, 1, 0, 3, 2],
        ).unwrap()
    }

    pub fn gen_tet_c() -> Manifold {
        Manifold::new(
            &vec![-2., -0.866025, 0.5, -2., -0., -1., -2., 0.866025, 0.5, 0., 0., 0.],
            &vec![1, 3, 0, 1, 0, 2, 2, 3, 1, 0, 3, 2],
        ).unwrap()
    }


    #[test]
    fn test_tet_sub_inclusion_case(){
        let expand = -1.;
        let mfd_p = gen_tet_a();
        let mfd_q = gen_tet_c();
        let mut p1q2 = vec![];
        let mut p2q1 = vec![];
        let (x12, v12) = intersect12(&mfd_p, &mfd_q, &mut p1q2, expand, true);
        let (x21, v21) = intersect12(&mfd_p, &mfd_q, &mut p2q1, expand, false);
        let w03 = winding03(&mfd_p, &mfd_q, expand, true);
        let w30 = winding03(&mfd_p, &mfd_q, expand, false);

        assert_eq!(w03, vec![0, 0, 0, 0]);
        assert_eq!(w30, vec![0, 0, 0, 1]);
        assert_eq!(x12.len(), 0);
        assert_eq!(v12.len(), 0);
        assert_eq!(x21, vec![-1, -1, -1]);
        let v21_ = vec![
            Vec3::new(-0.224009, 0., -0.112005),
            Vec3::new(-0.294367, 0.127465, 0.0735918),
            Vec3::new(-0.395087, -0.171077, 0.0987716),
        ];
        for i in 0..3 {
            assert!((v21[i] - v21_[i]).length() < 1e-6);
        }
        let op = OpType::Subtract;
        let b03 = Boolean03{ p1q2, p2q1, x12, x21, w03, w30, v12, v21 };
        let b45 = boolean45(&mfd_p, &mfd_q, &b03, &op);
    }

    #[test]
    fn test_tet_sub_penetration_case(){
        let expand = -1.;
        let mfd_p = gen_tet_a();
        let mfd_q = gen_tet_b();
        let mut p1q2 = vec![];
        let mut p2q1 = vec![];
        let (x12, v12) = intersect12(&mfd_p, &mfd_q, &mut p1q2, -1., true);
        let (x21, v21) = intersect12(&mfd_p, &mfd_q, &mut p2q1, -1., false);
        let w03 = winding03(&mfd_p, &mfd_q, expand, true);
        let w30 = winding03(&mfd_p, &mfd_q, expand, false);

        let v12_ = vec![
            Vec3::new(-0.763707, -0.763707, 0.440927),
            Vec3::new(-0.242656, 0.439609, 0.140098),
            Vec3::new(0.,0.,-0.5),
            Vec3::new(0.,0.,-0.5)
        ];
        let v21_ = vec![
            Vec3::new(0.302169,0.302169,0.174458),
            Vec3::new(0.439609,-0.242656,0.140098),
            Vec3::new(0.302169,0.302169,0.174458),
            Vec3::new(-0.763707,-0.763707,0.440927)
        ];

        assert_eq!(w03, vec![0, 0, 0, 0]);
        assert_eq!(w30, vec![0, 0, 0, 0]);
        assert_eq!(x12, vec![-1, 1, -1, 1]);
        assert_eq!(x21, vec![1, 1, -1, -1]);
        for i in 0..3 {
            assert!((v12[i] - v12_[i]).length() < 1e-6);
            assert!((v21[i] - v21_[i]).length() < 1e-6);
        }

        let op = OpType::Subtract;
        let b03 = Boolean03{ p1q2, p2q1, x12, x21, w03, w30, v12, v21 };
        let b45 = boolean45(&mfd_p, &mfd_q, &b03, &op);
    }
}

#[cfg(test)]
mod test_triangulation {
    use crate::triangulation::Pt;
    use crate::triangulation::ear_clip::EarClip;
    use crate::{Vec2, Vec3, Vec3u};

    #[test]
    fn test_ear_clip() {
        let polys = vec![
            vec![
                Pt {idx: 2120, pos: Vec2::new(0.048238, 0.680959)},
                Pt {idx: 2124, pos: Vec2::new(-0.0145625, -0.676874)},
                Pt {idx: 2123, pos: Vec2::new(0.0245192, -0.68213)},
                Pt {idx: 2119, pos: Vec2::new(0.0482562, -0.681659)},
            ],
            vec![
                Pt {idx: 2122, pos: Vec2::new(-0.068635, -0.673357)},
                Pt {idx: 2125, pos: Vec2::new(-0.0487738, -0.690778)},
                Pt {idx: 2121, pos: Vec2::new(-0.02279, -0.676339)},
            ],
        ];
        let res0 = vec![
            Vec3u::new(2123, 2119, 2120),
            Vec3u::new(2123, 2120, 2124),
            Vec3u::new(2125, 2121, 2122)
        ];

        let res1 = EarClip::new(&polys, 1e-12).triangulate();
        for i in 0..3 { assert_eq!(res0[i], res1[i]); }
    }
}

#[cfg(test)]
mod test_simplification {
    use crate::simplification::collapse::collapse_collinear_edges;
    use crate::Vec3;

    #[test]
    fn test_collapse() {
        use crate::common::{Half, Tref};

        let mut hs = vec![
            Half::new(0, 4, 5, usize::MAX),
            Half::new(4, 1, 42, usize::MAX),
            Half::new(1, 0, 9, usize::MAX),
            Half::new(0, 2, 26, usize::MAX),
            Half::new(2, 4, 50, usize::MAX),
            Half::new(4, 0, 0, usize::MAX),
            Half::new(0, 11, 11, usize::MAX),
            Half::new(11, 9, 57, usize::MAX),
            Half::new(9, 0, 21, usize::MAX),
            Half::new(0, 1, 2, usize::MAX),
            Half::new(1, 11, 14, usize::MAX),
            Half::new(11, 0, 6, usize::MAX),
            Half::new(1, 10, 17, usize::MAX),
            Half::new(10, 11, 68, usize::MAX),
            Half::new(11, 1, 10, usize::MAX),
            Half::new(1, 3, 35, usize::MAX),
            Half::new(3, 10, 20, usize::MAX),
            Half::new(10, 1, 12, usize::MAX),
            Half::new(3, 13, 31, usize::MAX),
            Half::new(13, 10, 66, usize::MAX),
            Half::new(10, 3, 16, usize::MAX),
            Half::new(0, 9, 8, usize::MAX),
            Half::new(9, 12, 59, usize::MAX),
            Half::new(12, 0, 24, usize::MAX),
            Half::new(0, 12, 23, usize::MAX),
            Half::new(12, 2, 27, usize::MAX),
            Half::new(2, 0, 3, usize::MAX),
            Half::new(2, 12, 25, usize::MAX),
            Half::new(12, 13, 64, usize::MAX),
            Half::new(13, 2, 30, usize::MAX),
            Half::new(2, 13, 29, usize::MAX),
            Half::new(13, 3, 18, usize::MAX),
            Half::new(3, 2, 36, usize::MAX),
            Half::new(1, 7, 47, usize::MAX),
            Half::new(7, 3, 37, usize::MAX),
            Half::new(3, 1, 15, usize::MAX),
            Half::new(2, 3, 32, usize::MAX),
            Half::new(3, 7, 34, usize::MAX),
            Half::new(7, 2, 51, usize::MAX),
            Half::new(4, 6, 49, usize::MAX),
            Half::new(6, 5, 54, usize::MAX),
            Half::new(5, 4, 43, usize::MAX),
            Half::new(1, 4, 1, usize::MAX),
            Half::new(4, 5, 41, usize::MAX),
            Half::new(5, 1, 45, usize::MAX),
            Half::new(1, 5, 44, usize::MAX),
            Half::new(5, 7, 56, usize::MAX),
            Half::new(7, 1, 33, usize::MAX),
            Half::new(2, 6, 53, usize::MAX),
            Half::new(6, 4, 39, usize::MAX),
            Half::new(4, 2, 4, usize::MAX),
            Half::new(2, 7, 38, usize::MAX),
            Half::new(7, 6, 55, usize::MAX),
            Half::new(6, 2, 48, usize::MAX),
            Half::new(5, 6, 40, usize::MAX),
            Half::new(6, 7, 52, usize::MAX),
            Half::new(7, 5, 46, usize::MAX),
            Half::new(9, 11, 7, usize::MAX),
            Half::new(11, 12, 61, usize::MAX),
            Half::new(12, 9, 22, usize::MAX),
            Half::new(8, 12, 65, usize::MAX),
            Half::new(12, 11, 58, usize::MAX),
            Half::new(11, 8, 69, usize::MAX),
            Half::new(8, 13, 71, usize::MAX),
            Half::new(13, 12, 28, usize::MAX),
            Half::new(12, 8, 60, usize::MAX),
            Half::new(10, 13, 19, usize::MAX),
            Half::new(13, 11, 70, usize::MAX),
            Half::new(11, 10, 13, usize::MAX),
            Half::new(8, 11, 62, usize::MAX),
            Half::new(11, 13, 67, usize::MAX),
            Half::new(13, 8, 63, usize::MAX),
        ];

        let mut ps = vec![
            Vec3::new(-0.8, -1., -1.),
            Vec3::new(-0.8, -1., 1.),
            Vec3::new(-0.8, 1., -1.),
            Vec3::new(-0.8, 1., 1.),
            Vec3::new(1.2, -1., -1.),
            Vec3::new(1.2, -1., 1.),
            Vec3::new(1.2, 1., -1.),
            Vec3::new(1.2, 1., 1.),
            Vec3::new(0., 0., 0.),
            Vec3::new(-0.8, -0.14641, -0.14641),
            Vec3::new(-0.8, 0.2, 0.2),
            Vec3::new(-0.8, -0.34641, 0.2),
            Vec3::new(-0.8, 0., -0.4),
            Vec3::new(-0.8, 0.34641, 0.2),
        ];

        let mut ns = vec![
            Vec3::new(-0., -1., -0.),
            Vec3::new(-0., -0., -1.),
            Vec3::new(-1., -0., 0.),
            Vec3::new(-1., -0., 0.),
            Vec3::new(-1., -0., 0.),
            Vec3::new(-1., -0., 0.),
            Vec3::new(-1., -0., -0.),
            Vec3::new(-1., -0., -0.),
            Vec3::new(-1., -0., -0.),
            Vec3::new(-1., -0., -0.),
            Vec3::new(-1., -0., -0.),
            Vec3::new(0., 0., 1.),
            Vec3::new(0., 1., 0.),
            Vec3::new(1., 0., 0.),
            Vec3::new(0., -1., 0.),
            Vec3::new(-0., 0., 1.),
            Vec3::new(0., 0., -1.),
            Vec3::new(-0., 1., 0.),
            Vec3::new(1., 0., -0.),
            Vec3::new(-0.242536, 0.840168, 0.485071),
            Vec3::new(-0.242536, 0.840168, 0.485071),
            Vec3::new(-0.242536, -0.840168, 0.485071),
            Vec3::new(-0.242536, -3.22877e-20, -0.970143),
            Vec3::new(-0.242536, -3.22877e-20, -0.970143),
        ];

        let mut refs = vec![
            Tref{mid: 1, fid: 0, pid: 3},
            Tref{mid: 1, fid: 0, pid: 5},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 2},
            Tref{mid: 1, fid: 0, pid: 1},
            Tref{mid: 1, fid: 0, pid: 0},
            Tref{mid: 1, fid: 0, pid: 4},
            Tref{mid: 1, fid: 0, pid: 3},
            Tref{mid: 1, fid: 0, pid: 1},
            Tref{mid: 1, fid: 0, pid: 5},
            Tref{mid: 1, fid: 0, pid: 0},
            Tref{mid: 1, fid: 0, pid: 4},
            Tref{mid: 2, fid: 0, pid: 0},
            Tref{mid: 2, fid: 0, pid: 0},
            Tref{mid: 2, fid: 0, pid: 3},
            Tref{mid: 2, fid: 0, pid: 2},
            Tref{mid: 2, fid: 0, pid: 2},
        ];

        collapse_collinear_edges(
            &mut hs,
            &mut ps,
            &mut ns,
            &mut refs,
            9,
            false,
            1e-6
        );

        let hs_out = vec![
            (0, 4, 5),
            (4, 1, 42),
            (1, 0, 9),
            (0, 2, 26),
            (2, 4, 50),
            (4, 0, 0),
            (0, 11, 11),
            (11, 12, 61),
            (12, 0, 24),
            (0, 1, 2),
            (1, 11, 17),
            (11, 0, 6),
            (-1, -1, -1),
            (-1, -1, -1),
            (-1, -1, -1),
            (1, 3, 35),
            (3, 11, 20),
            (11, 1, 10),
            (3, 13, 31),
            (13, 11, 70),
            (11, 3, 16),
            (-1, -1, -1),
            (-1, -1, -1),
            (-1, -1, -1),
            (0, 12, 8),
            (12, 2, 27),
            (2, 0, 3),
            (2, 12, 25),
            (12, 13, 64),
            (13, 2, 30),
            (2, 13, 29),
            (13, 3, 18),
            (3, 2, 36),
            (1, 7, 47),
            (7, 3, 37),
            (3, 1, 15),
            (2, 3, 32),
            (3, 7, 34),
            (7, 2, 51),
            (4, 6, 49),
            (6, 5, 54),
            (5, 4, 43),
            (1, 4, 1),
            (4, 5, 41),
            (5, 1, 45),
            (1, 5, 44),
            (5, 7, 56),
            (7, 1, 33),
            (2, 6, 53),
            (6, 4, 39),
            (4, 2, 4),
            (2, 7, 38),
            (7, 6, 55),
            (6, 2, 48),
            (5, 6, 40),
            (6, 7, 52),
            (7, 5, 46),
            (-1, -1, -1),
            (-1, -1, -1),
            (-1, -1, -1),
            (8, 12, 65),
            (12, 11, 7),
            (11, 8, 69),
            (8, 13, 71),
            (13, 12, 28),
            (12, 8, 60),
            (-1, -1, -1),
            (-1, -1, -1),
            (-1, -1, -1),
            (8, 11, 62),
            (11, 13, 19),
            (13, 8, 63),
        ];

        let ps_out = vec![
            Vec3::new(-0.8, -1., -1.),
            Vec3::new(-0.8, -1., 1.),
            Vec3::new(-0.8, 1., -1.),
            Vec3::new(-0.8, 1., 1.),
            Vec3::new(1.2, -1., -1.),
            Vec3::new(1.2, -1., 1.),
            Vec3::new(1.2, 1., -1.),
            Vec3::new(1.2, 1., 1.),
            Vec3::new(0., 0., 0.),
            Vec3::NAN,
            Vec3::NAN,
            Vec3::new(-0.8, -0.34641, 0.2),
            Vec3::new(-0.8, 0., -0.4),
            Vec3::new(-0.8, 0.34641, 0.2),
        ];

        for i in 0..hs.len() {
            let h = &hs[i];
            let a = hs[i].tail;
            let b = hs[i].head;
            let c = hs[i].pair;
            let (d, e, f) = hs_out[i];
            if h.tail().is_none() { assert_eq!(d, -1); } else { assert_eq!(a as i32, d); }
            if h.head().is_none() { assert_eq!(e, -1); } else { assert_eq!(b as i32, e); }
            if h.pair().is_none() { assert_eq!(f, -1); } else { assert_eq!(c as i32, f); }
        }

        for i in 0..ps.len() {
            let p0 = &ps[i];
            let p1 = &ps_out[i];
            if !p0.x.is_nan() && !p0.y.is_nan() && p0.z.is_nan() { assert!((p0 - p1).length() < 1e-6); }
        }
    }
}

#[cfg(test)]
mod test_mesh_cleanup {
    use crate::prelude::Manifold;

    #[test]
    fn test_dedup_verts() {
        let pos = vec![
           -0.866025, -1., 0.5, // duplicated
            0., -1., -1.,
            0.866025, -1., 0.5, // duplicated 2
           -0.866025, -1., 0.5, // duplicated
            0., 1., 0.,
            0.866025, -1., 0.5, // duplicated 2
        ];
        let idx = vec![
            0, 4, 1,
            0, 3, 1, // collapsed
            0, 3, 2, // collapsed
            0, 3, 4, // collapsed
            1, 2, 0,
            1, 4, 2,
            2, 4, 0,
            2, 5, 0, // collapsed
            5, 2, 1, // collapsed
        ];
        let mfd = Manifold::new(&pos, &idx).unwrap();

        assert_eq!(mfd.nv, 4);
        assert_eq!(mfd.nf, 4);
    }

    #[test]
    fn test_edge_shared_by_four_faces_is_an_error() {
        // two tetrahedra sharing the edge 0-1: that edge has four faces
        let pos = vec![
            0., 0., 0.,
            1., 0., 0.,
            0., 1., 0.,
            0., 0., 1.,
            0., -1., 0.,
            0., 0., -1.,
        ];
        let idx = vec![
            0, 2, 1,
            0, 1, 3,
            0, 3, 2,
            1, 2, 3,
            0, 4, 1,
            0, 1, 5,
            0, 5, 4,
            1, 4, 5,
        ];
        assert!(Manifold::new(&pos, &idx).is_err());
    }
}

// Degenerate configurations of two icospheres: identical, offset or rotated by a few ulp,
// tangent, nested, and differing in resolution. These exercise coincident faces and edges
// on curved surfaces, where |x12| == 2 entries can appear (see boolean45::duplicate_verts).
// Each case checks the volume identities  Add + Int == Va + Vb  and  Sub + Int == Va.
#[cfg(test)]
mod test_sphere_degenerate {
    use std::collections::HashMap;
    use crate::prelude::*;
    use crate::Vec3;

    fn icosphere(sub: usize, r: f64, c: [f64; 3], rot: f64) -> Manifold {
        let t = (1.0 + 5f64.sqrt()) / 2.0;
        let mut v: Vec<[f64; 3]> = vec![
            [-1., t, 0.], [1., t, 0.], [-1., -t, 0.], [1., -t, 0.],
            [0., -1., t], [0., 1., t], [0., -1., -t], [0., 1., -t],
            [t, 0., -1.], [t, 0., 1.], [-t, 0., -1.], [-t, 0., 1.],
        ];
        let mut f: Vec<[usize; 3]> = vec![
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
        ];
        for _ in 0..sub {
            let mut nf = vec![];
            let mut cache = HashMap::new();
            let mut mid = |a: usize, b: usize, v: &mut Vec<[f64; 3]>| -> usize {
                *cache.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    v.push([(v[a][0] + v[b][0]) / 2., (v[a][1] + v[b][1]) / 2., (v[a][2] + v[b][2]) / 2.]);
                    v.len() - 1
                })
            };
            for t in &f {
                let (a, b, c) = (t[0], t[1], t[2]);
                let ab = mid(a, b, &mut v);
                let bc = mid(b, c, &mut v);
                let ca = mid(c, a, &mut v);
                nf.extend_from_slice(&[[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
            }
            f = nf;
        }
        let (cs, sn) = (rot.cos(), rot.sin());
        let pos: Vec<f64> = v.iter().flat_map(|p| {
            let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let q = [p[0] / l * r, p[1] / l * r, p[2] / l * r];
            [cs * q[0] - sn * q[1] + c[0], sn * q[0] + cs * q[1] + c[1], q[2] + c[2]]
        }).collect();
        let idx: Vec<usize> = f.iter().flat_map(|t| t.iter().copied()).collect();
        Manifold::new(&pos, &idx).unwrap()
    }

    fn volume(m: &Manifold) -> f64 {
        m.hs.chunks(3).map(|t| {
            let a: Vec3 = m.ps[t[0].tail];
            let b: Vec3 = m.ps[t[1].tail];
            let c: Vec3 = m.ps[t[2].tail];
            a.dot(b.cross(c)) as f64
        }).sum::<f64>() / 6.
    }

    // An empty result is reported as Err("empty pos matrix"); treat it as zero volume.
    fn boolean_volume(a: &Manifold, b: &Manifold, op: OpType) -> f64 {
        match compute_boolean(a, b, op) {
            Ok(m) => volume(&m),
            Err(e) if e.contains("empty") => 0.,
            Err(e) => panic!("boolean failed: {e}"),
        }
    }

    fn check_identities(a: &Manifold, b: &Manifold) {
        let tol = if cfg!(feature = "f32") { 1e-4 } else { 1e-9 };
        let va = volume(a);
        let vb = volume(b);
        let add = boolean_volume(a, b, OpType::Add);
        let sub = boolean_volume(a, b, OpType::Subtract);
        let int = boolean_volume(a, b, OpType::Intersect);
        let e1 = (add + int - va - vb).abs() / (va + vb);
        let e2 = (sub + int - va).abs() / (va + vb);
        assert!(e1 < tol, "Add + Int != Va + Vb: add={add} int={int} va={va} vb={vb} rel_err={e1:e}");
        assert!(e2 < tol, "Sub + Int != Va: sub={sub} int={int} va={va} rel_err={e2:e}");
        assert!(int >= -tol && int <= va.min(vb) + tol, "Int out of range: {int}");
    }

    #[test]
    fn test_sphere_identical() {
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        check_identities(&s, &s);
        assert!(boolean_volume(&s, &s, OpType::Subtract).abs() < 1e-9);
    }

    #[test]
    fn test_sphere_shifted_by_ulps() {
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        for d in [1e-15, 1e-12, 1e-9, 1e-6, 1e-3] {
            check_identities(&s, &icosphere(3, 1., [d, 0., 0.], 0.));
        }
    }

    #[test]
    fn test_sphere_rotated_slightly() {
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        for d in [1e-12, 1e-9, 1e-6, 1e-3] {
            check_identities(&s, &icosphere(3, 1., [0., 0., 0.], d));
        }
    }

    // Known limitation: with a rotation of a few ulp the edges of the two spheres are
    // nearly coincident and the interpolation parameter in kernel01::intersect is
    // ill-conditioned, so intersection verts land at arbitrary points along the edges
    // and the two sheets of the result no longer enclose ~zero volume. Fixing this
    // needs adaptive-precision evaluation of the signed distances (Shewchuk-style).
    #[test]
    #[ignore]
    fn test_sphere_rotated_by_ulps() {
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        check_identities(&s, &icosphere(3, 1., [0., 0., 0.], 1e-15));
    }

    #[test]
    fn test_sphere_tangent_at_vertex() {
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        check_identities(&s, &icosphere(3, 1., [2., 0., 0.], 0.));
    }

    #[test]
    fn test_sphere_nested() {
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        check_identities(&s, &icosphere(3, 0.5, [0., 0., 0.], 0.));
    }

    #[test]
    fn test_sphere_vs_coarser_sphere() {
        // Shares the 12 icosahedron verts; every other vert of the coarser mesh lies inside the finer one.
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        check_identities(&s, &icosphere(2, 1., [0., 0., 0.], 0.));
    }

    #[test]
    fn test_sphere_rotated_36deg() {
        let s = icosphere(3, 1., [0., 0., 0.], 0.);
        check_identities(&s, &icosphere(3, 1., [0., 0., 0.], std::f64::consts::PI / 5.));
    }
}

#[cfg(test)]
mod test_values {
    use crate::prelude::*;
    use crate::{Real, Vec2, Vec3};
    use crate::values::get_barycentric_ratio;

    fn tol() -> f64 { if cfg!(feature = "f32") { 1e-4 } else { 1e-9 } }

    fn expected_uv(k: usize, s: f64, offset: f64, pos: [f64; 3]) -> [f64; 2] {
        let seam = if s > 0. { 0. } else { 5. };
        [pos[(k + 1) % 3] + offset + 10. * k as f64 + seam, pos[(k + 2) % 3] + offset]
    }

    fn uv_cube_buf(min: [f64; 3], size: f64, offset: f64) -> (Vec<f64>, Vec<usize>, Vec<Vec2>) {
        let (mut pos, mut idx, mut vals) = (vec![], vec![], vec![]);
        for k in 0..3 { for s in [1., -1.] {
            let (a, b) = ((k + 1) % 3, (k + 2) % 3);
            let base = pos.len() / 3;
            for (u, v) in [(0., 0.), (1., 0.), (1., 1.), (0., 1.)] {
                let mut p = min;
                if s > 0. { p[k] += size; }
                p[a] += u * size;
                p[b] += v * size;
                pos.extend_from_slice(&p);
                let uv = expected_uv(k, s, offset, p);
                vals.push(Vec2::new(uv[0] as Real, uv[1] as Real));
            }
            if s > 0. { idx.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]); }
            else      { idx.extend_from_slice(&[base, base + 2, base + 1, base, base + 3, base + 2]); }
        }}
        let vals = idx.iter().map(|&i| vals[i]).collect();
        (pos, idx, vals)
    }

    fn uv_cube(min: [f64; 3], size: f64, offset: f64) -> Manifold<Vec2> {
        let (pos, idx, vals) = uv_cube_buf(min, size, offset);
        Manifold::new_with_values(&pos, &idx, &vals).unwrap()
    }

    fn volume<S>(m: &Manifold<S>) -> f64 {
        m.hs.chunks(3).map(|t| {
            let (a, b, c) = (m.ps[t[0].tail], m.ps[t[1].tail], m.ps[t[2].tail]);
            a.dot(b.cross(c)) as f64
        }).sum::<f64>() / 6.
    }

    struct Cube { min: [f64; 3], size: f64, offset: Option<f64> }
    fn cube(min: [f64; 3], size: f64, offset: f64) -> Cube { Cube { min, size, offset: Some(offset) } }

    fn check_corner_uvs(m: &Manifold<Vec2>, cubes: &[Cube], tol: f64) {
        assert_eq!(m.vals.len(), m.hs.len());
        let corner = |i: usize| {
            let (p, uv) = (m.ps[m.hs[i].tail], m.vals[i]);
            ([p.x as f64, p.y as f64, p.z as f64], [uv.x as f64, uv.y as f64])
        };
        for t in 0..m.nf {
            let c = [corner(3 * t), corner(3 * t + 1), corner(3 * t + 2)];
            let k = (0..3)
                .find(|&k| (c[0].0[k] - c[1].0[k]).abs() < tol && (c[0].0[k] - c[2].0[k]).abs() < tol)
                .unwrap_or_else(|| panic!("tri {t} is not axis aligned: {c:?}"));
            let plane = c[0].0[k];
            let cen = [0, 1, 2].map(|a| (c[0].0[a] + c[1].0[a] + c[2].0[a]) / 3.);
            let mut cands = vec![];
            for (ci, cb) in cubes.iter().enumerate() {
                for s in [1., -1.] {
                    let f = cb.min[k] + if s > 0. { cb.size } else { 0. };
                    let on = (f - plane).abs() < tol;
                    let inside = (0..3).filter(|&a| a != k)
                        .all(|a| cen[a] > cb.min[a] - tol && cen[a] < cb.min[a] + cb.size + tol);
                    if on && inside { cands.push((ci, s)); }
                }
            }
            assert!(!cands.is_empty(), "tri {t} on axis {k} = {plane} with centroid {cen:?} lies on no cube face");
            for (i, (p, uv)) in c.iter().enumerate() {
                let ok = cands.iter().any(|&(ci, s)| {
                    let e = match cubes[ci].offset { Some(o) => expected_uv(k, s, o, *p), None => [0., 0.] };
                    (e[0] - uv[0]).abs() < tol && (e[1] - uv[1]).abs() < tol
                });
                assert!(ok, "tri {t} corner {i} at {p:?}: uv {uv:?} matches none of the faces {cands:?} (axis {k} = {plane})");
            }
        }
    }

    #[test]
    fn test_values_input() {
        let m = uv_cube([0., 0., 0.], 1., 0.);
        assert_eq!((m.nv, m.nf), (8, 12));
        assert_eq!(m.vals.len(), 36);
        assert!((volume(&m) - 1.).abs() < tol());
        let mut per_vert = vec![std::collections::BTreeSet::new(); m.nv];
        for h in &m.hs { per_vert[h.tail].insert(h.prop); }
        assert!(per_vert.iter().all(|s| s.len() == 3));
        check_corner_uvs(&m, &[cube([0., 0., 0.], 1., 0.)], tol());
    }

    #[test]
    fn test_values_none() {
        let (pos, idx, _) = uv_cube_buf([0., 0., 0.], 1., 0.);
        let m = Manifold::<Vec2>::new_with_values(&pos, &idx, &[]).unwrap();
        assert!(m.vals.is_empty());
        assert!(m.hs.iter().all(|h| h.prop == h.tail));
        let n = uv_cube([0., 0., 0.], 1., 0.);
        assert_eq!((n.nv, n.nf, n.nh), (m.nv, m.nf, m.nh));
        assert_eq!(n.ps, m.ps);
        for (a, b) in n.hs.iter().zip(m.hs.iter()) {
            assert_eq!((a.tail, a.head, a.pair), (b.tail, b.head, b.pair));
        }
    }

    #[test]
    fn test_values_scalar() {
        let f = |p: &[f64]| (p[0] + 2. * p[1] - 3. * p[2]) as Real;
        let (pos, idx, _) = uv_cube_buf([0., 0., 0.], 1., 0.);
        let p = Manifold::new_with_values(&pos, &idx, &idx.iter().map(|&i| f(&pos[3 * i..3 * i + 3])).collect::<Vec<_>>()).unwrap();
        let (pos, idx, _) = uv_cube_buf([0.3, 0.4, 0.5], 1., 0.);
        let q = Manifold::new_with_values(&pos, &idx, &idx.iter().map(|&i| f(&pos[3 * i..3 * i + 3])).collect::<Vec<_>>()).unwrap();
        for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
            let r = compute_boolean(&p, &q, op).unwrap();
            for (h, &v) in r.hs.iter().zip(&r.vals) {
                let x = r.ps[h.tail];
                let e = f(&[x.x as f64, x.y as f64, x.z as f64]);
                assert!(((v - e) as f64).abs() < tol(), "{v} != {e}");
            }
        }
    }

    #[test]
    fn test_values_disjoint_union() {
        let p = uv_cube([0., 0., 0.], 1., 0.);
        let q = uv_cube([3., 0., 0.], 1., 40.);
        let r = compute_boolean(&p, &q, OpType::Add).unwrap();
        assert!(!r.vals.is_empty());
        assert!((volume(&r) - 2.).abs() < tol());
        check_corner_uvs(&r, &[cube([0., 0., 0.], 1., 0.), cube([3., 0., 0.], 1., 40.)], tol());
    }

    #[test]
    fn test_values_overlapping() {
        let p = uv_cube([0., 0., 0.], 1., 0.);
        let q = uv_cube([0.5, 0.5, 0.5], 1., 40.);
        let cubes = [cube([0., 0., 0.], 1., 0.), cube([0.5, 0.5, 0.5], 1., 40.)];
        for (op, vol) in [(OpType::Add, 2. - 0.125), (OpType::Subtract, 1. - 0.125), (OpType::Intersect, 0.125)] {
            let r = compute_boolean(&p, &q, op).unwrap();
            assert!(!r.vals.is_empty());
            assert!((volume(&r) - vol).abs() < tol(), "volume {} != {vol}", volume(&r));
            check_corner_uvs(&r, &cubes, tol());
        }
    }

    #[test]
    fn test_values_coplanar_simplification() {
        for (min, add, sub) in [([0.5, 0., 0.], 1.5, 0.5), ([1., 0., 0.], 2., 1.), ([0.5, 0.5, 0.], 1.75, 0.75)] {
            let p = uv_cube([0., 0., 0.], 1., 0.);
            let q = uv_cube(min, 1., 40.);
            let cubes = [cube([0., 0., 0.], 1., 0.), cube(min, 1., 40.)];
            for (op, vol) in [(OpType::Add, add), (OpType::Subtract, sub)] {
                let r = compute_boolean(&p, &q, op).unwrap();
                assert!(!r.vals.is_empty());
                assert!((volume(&r) - vol).abs() < tol(), "min = {min:?}: volume {} != {vol}", volume(&r));
                check_corner_uvs(&r, &cubes, tol());
            }
        }
    }

    #[test]
    fn test_values_mismatch() {
        let p = uv_cube([0., 0., 0.], 1., 0.);
        let (pos, idx, _) = uv_cube_buf([0.5, 0.5, 0.5], 1., 0.);
        let q = Manifold::new_with_values(&pos, &idx, &[]).unwrap();
        let cubes = [cube([0., 0., 0.], 1., 0.), Cube { min: [0.5, 0.5, 0.5], size: 1., offset: None }];
        for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
            let r = compute_boolean(&p, &q, op).unwrap();
            assert!(!r.vals.is_empty());
            check_corner_uvs(&r, &cubes, tol());
        }
    }

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

#[cfg(test)]
mod test_input_validation {
    use crate::prelude::Manifold;

    // A unit tetrahedron: valid input, so each test changes one thing.
    const POS: [f64; 12] = [0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1.];
    const IDX: [usize; 12] = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];

    #[test]
    fn test_valid_tet_is_accepted() {
        assert!(Manifold::new(&POS, &IDX).is_ok());
    }

    #[test]
    fn test_index_past_end_of_pos_is_an_error() {
        let mut idx = IDX;
        idx[0] = 99;
        assert!(Manifold::new(&POS, &idx).is_err());
    }

    #[test]
    fn test_index_past_last_vertex_is_an_error() {
        // 4 vertices, but `pos.len()` is 12: index 5 used to be read as vertex 0.
        let mut idx = IDX;
        idx[0] = 5;
        assert!(Manifold::new(&POS, &idx).is_err());
    }

    #[test]
    fn test_non_finite_coordinates_are_errors() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut pos = POS;
            pos[3] = bad;
            assert!(Manifold::new(&pos, &IDX).is_err(), "{bad} was accepted");
        }
    }
}
