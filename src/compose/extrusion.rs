use std::f64::consts::PI;

use geo::{BoundingRect, Coord, LineString, MultiPolygon, Polygon};
use thiserror::Error;

use crate::{
    common::{Affine3, Vec3u},
    manifold::ManifoldError,
    prelude::Manifold,
    triangulation::{ear_clip::EarClip, Pt},
    Mat2, Mat3, Real, Vec2, Vec3, K_PRECISION,
};

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

/// Controls how faces are created on an extruded/revolved manifold.
enum FaceMode {
    /// Close the face
    Close,

    /// No face, loop the structure back to its first layer
    Loop,
}

fn raw_extrude_impl<'s, P>(
    polygon_iter_builder: impl Fn() -> P,
    divisions: usize,
    face_mode: FaceMode,
    affine: impl Fn(Real) -> Affine3,
) -> Result<Manifold, ManifoldError>
where
    P: IntoIterator<Item = &'s Polygon<Real>>,
{
    fn points(polygon: &Polygon<Real>) -> impl Iterator<Item = Vec3> {
        polygon
            .coords()
            .map(|coord| Vec3::new(coord.x, coord.y, 0.0))
    }

    let face_indicies = if matches!(face_mode, FaceMode::Close) {
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

        Some(EarClip::new(&polygons, K_PRECISION).triangulate())
    } else {
        // If we don't need to close the faces, then we don't need to calculate the face indicies.
        None
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

    if let Some(face_indicies) = face_indicies.as_ref() {
        // Insert bottom vertex references.
        for i in face_indicies.iter() {
            oft_ts.push(Vec3u::new(i.z, i.y, i.x));
        }
    }

    // Incert divisions. Note that the top of the shape counts as a division.
    for layer in 0..divisions {
        let alpha = (layer + 1) as Real / divisions as Real;
        let affine = affine(alpha);

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

    if let Some(face_indicies) = face_indicies.as_ref() {
        // Insert top vertex references.
        // We do not need to insert their verticies because they were provided by the final layer of
        // the divisions loop.
        for i in face_indicies.iter() {
            oft_ts.push(Vec3u::new(
                i.x + points_per_division * divisions,
                i.y + points_per_division * divisions,
                i.z + points_per_division * divisions,
            ));
        }
    } else {
        // Loop the final layer back to the first layer.
        let mut polygon_point_offset = 0;
        for polygon in polygon_iter_builder() {
            let ending_offset = points_per_division * divisions + polygon_point_offset;
            let starting_offset = polygon_point_offset;
            let points_in_polygon = polygon.num_coords();
            for (vertex_index, _position) in points(polygon).enumerate() {
                // Corners of a quardrangle making up a a side of the extruded shape.
                // k--l
                // |  |
                // i--j
                let k = vertex_index + starting_offset;
                let l = (vertex_index + 1) % points_in_polygon + starting_offset;
                let i = vertex_index + ending_offset;
                let j = (vertex_index + 1) % points_in_polygon + ending_offset;

                oft_ts.push(Vec3u::new(i, j, k));
                oft_ts.push(Vec3u::new(k, j, l));
            }

            polygon_point_offset += points_in_polygon;
        }
    }

    Manifold::new_impl(oft_ps, oft_ts, None, None)
}

#[derive(Debug, Error)]
pub enum ExtrusionError {
    #[error("Extrusion height must be greater than zero")]
    InvalidHeight,

    #[error("Height of extrusion top must be greater than or equel to zero")]
    InvalidScale,

    #[error("Error buiding manifold: {0}")]
    Manifold(#[from] ManifoldError),
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

    let manifold = raw_extrude_impl(polygon_iter_builder, divisions, FaceMode::Close, |alpha| {
        let phi = alpha * twist_radians;
        let scale = Vec2::splat(1.0).lerp(scale_top, alpha);
        let translation = Vec3::new(0.0, 0.0, alpha * height);
        let matrix2 = Mat2::from_scale_angle(scale, phi);
        let matrix3 = Mat3::from_mat2(matrix2);
        Affine3 {
            matrix3,
            translation,
        }
    })?;

    Ok(manifold)
}

#[derive(Debug, Error)]
pub enum RevolveError {
    #[error("Revolution angle must be greater than zero")]
    InvalidAngle,

    #[error("Geometry must not be present on the left side of the Y axis")]
    LeftOfYAxis,

    #[error("Error buiding manifold: {0}")]
    Manifold(#[from] ManifoldError),
}

fn revolve_impl<'s, P>(
    polygon_iter_builder: impl Fn() -> P,
    divisions: usize,
    angle_radians: Real,
) -> Result<Manifold, RevolveError>
where
    P: IntoIterator<Item = &'s Polygon<Real>>,
{
    if angle_radians <= 0.0 {
        return Err(RevolveError::InvalidAngle);
    }

    // Checks if any geometry has a point left of the Y axis
    if polygon_iter_builder().into_iter().any(|polygon| {
        polygon
            .bounding_rect()
            .is_some_and(|rect| rect.min().x < 0.0)
    }) {
        return Err(RevolveError::LeftOfYAxis);
    }

    const MAX_ANGLE: Real = PI as Real * 2.0;

    // Cap the angle at 2Pi.
    let angle_radians = MAX_ANGLE.min(angle_radians);

    let face_mode = if angle_radians < MAX_ANGLE {
        FaceMode::Close
    } else {
        FaceMode::Loop
    };

    let manifold = raw_extrude_impl(polygon_iter_builder, divisions, face_mode, |alpha| {
        let translation = Vec3::new(0.0, 0.0, 0.0);

        #[allow(clippy::useless_conversion)]
        let matrix3 = Mat3::from_axis_angle(Vec3::Y.into(), -angle_radians * alpha);
        Affine3 {
            matrix3,
            translation,
        }
    })?;

    Ok(manifold)
}

pub trait ExtrudePoly {
    fn extrude(
        &self,
        height: Real,
        divisions: usize,
        twist_radians: Real,
        scale_top: Vec2,
    ) -> Result<Manifold, ExtrusionError>;

    fn revolve(&self, divisions: usize, angle_radians: Real) -> Result<Manifold, RevolveError>;
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

    fn revolve(&self, divisions: usize, angle_radians: Real) -> Result<Manifold, RevolveError> {
        revolve_impl(|| [self], divisions, angle_radians)
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

    fn revolve(&self, divisions: usize, angle_radians: Real) -> Result<Manifold, RevolveError> {
        revolve_impl(|| self.iter(), divisions, angle_radians)
    }
}
