use geo::{Coord, LineString, MultiPolygon, Polygon};
use thiserror::Error;

use crate::{
    common::{Affine3, Vec3u},
    manifold::ManifoldError,
    prelude::Manifold,
    triangulation::{ear_clip::EarClip, Pt},
    Mat2, Mat3, Real, Vec2, Vec3, K_PRECISION,
};

#[derive(Debug, Error)]
pub enum ExtrusionError {
    #[error("Extrusion height must be greater than zero")]
    InvalidHeight,

    #[error("Height of extrusion top must be greater than or equel to zero")]
    InvalidScale,

    #[error("Error buiding manifold: {0}")]
    Manifold(#[from] ManifoldError),
}

pub trait ExtrudePoly {
    fn extrude(
        &self,
        height: Real,
        divisions: usize,
        twist_radians: Real,
        scale_top: Vec2,
    ) -> Result<Manifold, ExtrusionError>;
}

trait IterStrings {
    fn strings(&self) -> impl Iterator<Item = &LineString<Real>>;
    fn coords(&self) -> impl Iterator<Item = &Coord<Real>>;
    fn num_coords(&self) -> usize;
}

impl IterStrings for Polygon<Real> {
    fn strings(&self) -> impl Iterator<Item = &LineString<Real>> {
        [self.exterior()].into_iter().chain(self.interiors())
    }

    fn coords(&self) -> impl Iterator<Item = &Coord<Real>> {
        self.strings().flat_map(|string| string.coords())
    }

    fn num_coords(&self) -> usize {
        // TODO instead of counting each coor individually, sum up the length of all the strings.
        self.coords().count()
    }
}

fn extrude_impl<'s, P>(
    polygon_iter_builder: impl Fn() -> P,
    height: Real,
    divisions: usize,
    twist_radians: Real,
    scale_top: Vec2,
) -> Result<Manifold, ExtrusionError>
where
    P: IntoIterator<Item = &'s Polygon<Real>>,
{
    if height <= 0.0 {
        return Err(ExtrusionError::InvalidHeight);
    }

    if scale_top.x < 0.0 || scale_top.y < 0.0 {
        return Err(ExtrusionError::InvalidScale);
    }

    fn points(polygon: &Polygon<Real>) -> impl Iterator<Item = Vec3> {
        polygon
            .coords()
            .map(|coord| Vec3::new(coord.x, coord.y, 0.0))
    }

    let idcs = {
        let mut i = 0;

        // TODO I dislike collecting the polygons into a throw-away vec like this, but I want to
        // avoid changing the core library for now.
        let polygons: Vec<Vec<Pt>> = polygon_iter_builder()
            .into_iter()
            .map(|polygon| {
                polygon
                    .coords()
                    .map(|c| {
                        let pt = Pt {
                            pos: Vec2::new(c.x, c.y),
                            idx: i,
                        };
                        i += 1;
                        pt
                    })
                    .collect::<Vec<Pt>>()
            })
            .collect();

        EarClip::new(&polygons, K_PRECISION).triangulate()
    };

    let mut oft_ps = vec![];
    let mut oft_ts = vec![];
    let points_per_division: usize = polygon_iter_builder()
        .into_iter()
        .map(|polygon| polygon.num_coords())
        .sum();

    // Insert bottom verticies.
    for p in polygon_iter_builder().into_iter().flat_map(points) {
        oft_ps.push(p);
    }

    // Insert bottom vertex references.
    for i in idcs.iter() {
        oft_ts.push(Vec3u::new(i.z, i.y, i.x));
    }

    // Incert divisions. Note that the top of the shape counts as a division.
    for layer in 0..divisions {
        let alpha = (layer + 1) as Real / divisions as Real;

        // TODO we should let an external function calculate the affine.
        let phi = alpha * twist_radians;
        let scale = Vec2::splat(1.0).lerp(scale_top, alpha);
        let translation = Vec3::new(0.0, 0.0, alpha * height);
        let matrix2 = Mat2::from_scale_angle(scale, phi);
        let matrix3 = Mat3::from_mat2(matrix2);
        let affine = Affine3 {
            matrix3,
            translation,
        };

        // Insert the next division's verticies
        let mut polygon_point_offset = 0;
        for polygon in polygon_iter_builder() {
            let base_offset = layer * points_per_division + polygon_point_offset;
            let points_in_polygon = polygon.num_coords();
            for (vertex_index, position) in points(polygon).enumerate() {
                // Conversion is necessary for 32bit support.
                #[allow(clippy::useless_conversion)]
                oft_ps.push(affine.transform_point3(position.into()).into());

                // Corners of a quardrangle making up a a side of the extruded shape.
                // k--l
                // |  |
                // i--j
                let i = base_offset + vertex_index;
                let j = base_offset + (vertex_index + 1) % points_in_polygon;
                let k = i + points_per_division;
                let l = j + points_per_division;

                oft_ts.push(Vec3u::new(i, j, k));
                oft_ts.push(Vec3u::new(k, j, l));
            }

            polygon_point_offset += points_in_polygon;
        }
    }

    // Insert top vertex references.
    // We do not need to insert their verticies because they were provided by the final layer of
    // the divisions loop.
    for i in idcs.iter() {
        oft_ts.push(Vec3u::new(
            i.x + points_per_division * divisions,
            i.y + points_per_division * divisions,
            i.z + points_per_division * divisions,
        ));
    }

    Ok(Manifold::new_impl(oft_ps, oft_ts, None, None)?)
}

impl ExtrudePoly for Polygon<Real> {
    fn extrude(
        &self,
        height: Real,
        divisions: usize,
        twist_radians: Real,
        scale_top: Vec2,
    ) -> Result<Manifold, ExtrusionError> {
        extrude_impl(|| [self], height, divisions, twist_radians, scale_top)
    }
}

impl ExtrudePoly for MultiPolygon<Real> {
    fn extrude(
        &self,
        height: Real,
        divisions: usize,
        twist_radians: Real,
        scale_top: Vec2,
    ) -> Result<Manifold, ExtrusionError> {
        extrude_impl(|| self.iter(), height, divisions, twist_radians, scale_top)
    }
}
