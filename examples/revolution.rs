//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use std::f64::consts::PI;

use bevy::asset::RenderAssetUsages;
use bevy::color::palettes::css::*;
use bevy::pbr::wireframe::{Wireframe, WireframeColor, WireframePlugin};
use bevy::prelude::*;
use bevy::render::mesh::{PrimitiveTopology, VertexAttributeValues};
use bevy_panorbit_camera::{PanOrbitCamera, PanOrbitCameraPlugin};
use boolmesh::{prelude::*, Real};
use geo::{BooleanOps, Coord, MultiPolygon, Rect, Translate};

#[derive(Default, Reflect, GizmoConfigGroup)]
struct MyRoundGizmos {}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(WireframePlugin::default())
        .add_plugins(PanOrbitCameraPlugin)
        .init_gizmo_group::<MyRoundGizmos>()
        .add_systems(Startup, setup)
        .run();
}

fn setup(
    mut cmds: Commands,
    mut mats: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    fn square_with_a_bite_out_of_it() -> MultiPolygon<Real> {
        let square = Rect::new(Coord { x: 0.0, y: -0.5 }, Coord { x: 1.0, y: 0.5 }).to_polygon();
        let bite = Rect::new(Coord { x: 0.25, y: -0.25 }, Coord { x: 0.75, y: 0.5 }).to_polygon();
        square.boolean_op(&bite, geo::OpType::Difference)
    }

    fn two_polygons() -> MultiPolygon<Real> {
        let square = Rect::new(Coord { x: -0.5, y: -0.5 }, Coord { x: 0.5, y: 0.5 }).to_polygon();
        let square2 =
            Rect::new(Coord { x: -0.25, y: -0.25 }, Coord { x: 0.25, y: 0.25 }).to_polygon();
        let square2_hole = Rect::new(
            Coord {
                x: -0.124,
                y: -0.125,
            },
            Coord { x: 0.125, y: 0.125 },
        )
        .to_polygon();
        let square2 = square2.difference(&square2_hole).translate(0.0, 2.0);

        square
            .boolean_op(&square2, geo::OpType::Union)
            .translate(0.5, 0.0)
    }

    spawn_model(
        &mut cmds,
        &mut meshes,
        &mut mats,
        Vec3::new(-6.0, 0.0, 0.0),
        || square_with_a_bite_out_of_it().translate(0.5, 0.0),
        |polygon| polygon.revolve(5, PI as Real * 0.5).unwrap(),
    );

    spawn_model(
        &mut cmds,
        &mut meshes,
        &mut mats,
        Vec3::new(-3.0, 0.0, 0.0),
        square_with_a_bite_out_of_it,
        |polygon| polygon.revolve(15, PI as Real * 2.0).unwrap(),
    );

    spawn_model(
        &mut cmds,
        &mut meshes,
        &mut mats,
        Vec3::new(0.0, 0.0, 0.0),
        || square_with_a_bite_out_of_it().translate(0.5, 0.0),
        |polygon| polygon.revolve(15, PI as Real * 2.0).unwrap(),
    );

    spawn_model(
        &mut cmds,
        &mut meshes,
        &mut mats,
        Vec3::new(3.0, 0.0, 0.0),
        two_polygons,
        |polygon| polygon.revolve(15, PI as Real * 2.0).unwrap(),
    );

    spawn_model(
        &mut cmds,
        &mut meshes,
        &mut mats,
        Vec3::new(6.0, 0.0, 0.0),
        two_polygons,
        |polygon| polygon.revolve(5, PI as Real * 0.5).unwrap(),
    );

    cmds.spawn((PointLight::default(), Transform::from_xyz(2., 5., 2.)));
    cmds.spawn((
        Transform::from_translation(Vec3::new(0., 2., 3.)),
        PanOrbitCamera::default(),
    ));
}

fn spawn_model(
    cmds: &mut Commands,
    meshes: &mut ResMut<Assets<Mesh>>,
    mats: &mut ResMut<Assets<StandardMaterial>>,
    translation: Vec3,
    build_polygon: impl FnOnce() -> MultiPolygon<Real>,
    build_model: impl FnOnce(MultiPolygon<Real>) -> Manifold,
) {
    let polygon = build_polygon();

    let polygon_meshes = geo_bevy::multi_polygon_to_mesh(&polygon).unwrap();

    for mesh in polygon_meshes {
        let mut mesh = mesh.mesh;
        if let VertexAttributeValues::Float32x3(values) =
            mesh.attribute_mut(Mesh::ATTRIBUTE_NORMAL).unwrap()
        {
            values.iter_mut().for_each(|value| *value = [0.0, 0.0, 1.0]);
        }

        cmds.spawn((
            Mesh3d(meshes.add(mesh).clone()),
            MeshMaterial3d(mats.add(StandardMaterial { ..default() })),
            Transform::from_translation(translation + Vec3::new(0.0, 0.0, 2.0)),
            Wireframe,
            WireframeColor {
                color: BLACK.into(),
            },
        ));
    }

    let model = build_model(polygon);

    let mut m = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    let mut pos = vec![];
    let mut vns = vec![];
    for (fid, hs) in model.hs.chunks(3).enumerate() {
        let p0 = model.ps[hs[0].tail];
        let p1 = model.ps[hs[1].tail];
        let p2 = model.ps[hs[2].tail];
        let n = model.face_normals[fid];
        pos.push([p0.x as f32, p0.y as f32, p0.z as f32]);
        pos.push([p1.x as f32, p1.y as f32, p1.z as f32]);
        pos.push([p2.x as f32, p2.y as f32, p2.z as f32]);
        vns.push([n.x as f32, n.y as f32, n.z as f32]);
        vns.push([n.x as f32, n.y as f32, n.z as f32]);
        vns.push([n.x as f32, n.y as f32, n.z as f32]);
    }
    m.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    m.insert_attribute(Mesh::ATTRIBUTE_NORMAL, vns);

    cmds.spawn((
        Mesh3d(meshes.add(m).clone()),
        MeshMaterial3d(mats.add(StandardMaterial { ..default() })),
        Transform::from_translation(translation),
        Wireframe,
        WireframeColor {
            color: BLACK.into(),
        },
    ));
}
