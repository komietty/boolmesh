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

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::VecDeque;

use geo::{LineString, MultiPolygon, Polygon, Coord};
use thiserror::Error;

use crate::boolean03::boolean03;
use crate::boolean45::boolean45;
use crate::common::*;
use crate::manifold::bounds::Query;
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

#[derive(Debug, Error)]
pub enum ProjectionError {
    #[error("No polygons were produced by the operation")]
    NoPolygons,
}

/// Projects the manifold onto the XY plane. Rotate the manifold to project onto custom planes.
/// projection.
/// * manifold - Input manifold to project
pub fn compute_projection(manifold: &Manifold) -> Result<MultiPolygon, ProjectionError> {
    // TODO there should be a way to directly iterate triangles.
    let mut edge_ids: BTreeMap<usize, VecDeque<usize>> = BTreeMap::new();

    trait EdgeMap {
        fn next_starting_edge(&self) -> Option<usize>;
        fn next_edge(&mut self, manifold: &Manifold, current_edge_id: usize) -> Option<usize>;
    }

    impl EdgeMap for BTreeMap<usize, VecDeque<usize>> {
        fn next_starting_edge(&self) -> Option<usize> {
            let (_edge_id, queue) = self.first_key_value()?;
            let value = queue.back();
            value.copied()
        }
        
        fn next_edge(&mut self, manifold: &Manifold, current_edge_id: usize) -> Option<usize> {
            let current_key = manifold.hs[current_edge_id].head;
            let queue = self.get_mut(&current_key)?;
            let value = queue.pop_back();
            if queue.is_empty() {
                self.remove(&current_key);
            }

            value
        }
    }

    for (edge_id, edge) in manifold.hs.iter().enumerate() {
        // This filters our faces so that only faces that are connected to another face that is on
        // the opposite side of the manifold are included. This instantly gives us the edge
        // boundaries.
        if manifold.face_normals[manifold.hs[edge.pair].pair / 3].z <= 0.0
         && manifold.face_normals[edge.pair / 3].z > 0.0 {
            edge_ids.entry(edge.tail).or_default().push_front(edge_id);
        }
    }

    let mut polygons = Vec::new();
    while let Some(first_edge_id) = edge_ids.next_starting_edge() {
        let mut current_edge_id = first_edge_id;
        let mut line_string = Vec::new();

        loop {
            let point = manifold.ps[manifold.hs[current_edge_id].head];
            line_string.push(Coord { x: point.x, y: point.y });

            let next_edge_id =  edge_ids.next_edge(manifold, current_edge_id).expect("Non-manafold edge");
            
            if next_edge_id != first_edge_id {
                current_edge_id = next_edge_id;
            } else {
                // We've come back to our initial point.
                break;
            }
        }

        let mut line_string = LineString(line_string);
        line_string.close();
        let raw_polygon = Polygon::new(line_string, vec![]);

        polygons.push(raw_polygon);
    }

    let polygon = geo::unary_union(&polygons);

    if polygons.is_empty() {
        Err(ProjectionError::NoPolygons)
    } else {
        Ok(polygon)
    }
}

#[derive(Debug, Error)]
pub enum SliceError {
    #[error("No polygons were produced by the operation")]
    NoPolygons,
}

/// Slice a manifold into a 2D polygon
/// * manifold - Input manifold to slice
/// * height - z height to slice at
pub fn compute_slice(manifold: &Manifold, height: Real) -> Result<MultiPolygon, SliceError> {
    let mut bounding_box = manifold.bounding_box.clone();
    bounding_box.min.z = height;
    bounding_box.max.z = height;
    bounding_box.id = Some(0); // Collider will not report collisions without this.

    let mut triangle_ids = BTreeSet::new();

    manifold
        .collider
        .collision(&[Query::Bb(bounding_box)], &mut |_query_id, triangle_id| {
            let z_points = [0, 1, 2]
                .into_iter()
                .map(|j| manifold.ps[manifold.hs[3 * triangle_id + j].tail].z);

            // We have to account for NaN with these min/max functions, so we're going to just
            // filter out the NaNs.
            let min = z_points
                .clone()
                .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Greater));
            let max = z_points.max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less));

            // If the lowest point is below the height threashold, and the highest point is above,
            // then this triangle intersects with the height plane.
            if let (Some(min), Some(max)) = (min, max) && min <= height && max > height {
                triangle_ids.insert(triangle_id);
            }
        });

    // At this point, triangle_ids contains a list of triangles that intersect with the height
    // plane.
    fn next3(j: usize) -> usize {
      (j + 1) % 3
    }

    let mut polygons = Vec::new();

    while !triangle_ids.is_empty() {
        let start_triangle_id = *triangle_ids.first().ok_or(SliceError::NoPolygons)?;
        
        let mut vertex_index = 0;
        for j in [0, 1, 2] {
            if manifold.ps[manifold.hs[3 * start_triangle_id + j].tail].z > height &&
                manifold.ps[manifold.hs[3 * start_triangle_id + next3(j)].tail].z <= height {
              vertex_index = next3(j);
              break;
            }
        }

        let mut line_string = Vec::new();
        let mut current_triangle_id = start_triangle_id;
        loop {
            triangle_ids.remove(&current_triangle_id);

            if manifold.ps[manifold.hs[3 * current_triangle_id + vertex_index].head].z <= height {
              vertex_index = next3(vertex_index);
            }

            let up = &manifold.hs[3 * current_triangle_id + vertex_index];
            let below = manifold.ps[up.tail];
            let above = manifold.ps[up.head];
            let a = (height - below.z) / (above.z - below.z);
            let point = below.lerp(above, a);
            line_string.push(geo::Coord { x: point.x, y: point.y });

            let pair = up.pair;
            current_triangle_id = pair / 3;
            vertex_index = next3(pair % 3);

            if current_triangle_id == start_triangle_id {
                break;
            }
        }

        let mut line_string = LineString(line_string);
        line_string.close();
        polygons.push(Polygon::new(line_string, vec![]));
    }
   
    let polygon = geo::unary_union(&polygons);

    if polygons.is_empty() {
        Err(SliceError::NoPolygons)
    } else {
        Ok(polygon)
    }
}
