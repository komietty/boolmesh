//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use crate::common::{Half, Tref};
use crate::manifold::bounds::BBox;
use crate::manifold::collider::{morton_code, K_NO_CODE};
use crate::Vec3;

pub fn cleanup_unused_verts(
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
    *rs = hs.chunks(3).enumerate().filter_map(|(i, t)| t[0].pair().map(|_| rs[i].clone())).collect();
    *hs = hs.iter().filter(|h| h.pair().is_some()).cloned().collect();
}
