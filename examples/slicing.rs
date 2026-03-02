//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use std::f64::consts::PI;

use bevy::asset::RenderAssetUsages;
use bevy::color::palettes::css::*;
use bevy::pbr::wireframe::{Wireframe, WireframeColor, WireframePlugin};
use bevy::prelude::*;
use bevy::render::mesh::{PrimitiveTopology, VertexAttributeValues};
use bevy_panorbit_camera::{PanOrbitCamera, PanOrbitCameraPlugin};
use boolmesh::{compute_slice, prelude::*, Real};

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
    let model = "examples/models/double-torus.obj";
    let (gargoyle, _) = tobj::load_obj(
        model,
        &tobj::LoadOptions {
            ..Default::default()
        },
    )
    .expect("Failed to load obj file");

    let m = &gargoyle[0].mesh;

    let mut model = Manifold::new(
        &m.positions.iter().map(|&v| v as f64).collect::<Vec<_>>(),
        &m.indices.iter().map(|&v| v as usize).collect::<Vec<_>>(),
    )
    .unwrap();

    model.translate(-6.0, 0.0, 0.0);
    for _ in 0..5 {
        model.rotate(PI / 4.0, 0.0, 0.0);
        model.translate(2.0, 0.0, 0.0);

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
            Transform::default(),
            Wireframe,
            WireframeColor {
                color: BLACK.into(),
            },
        ));

        for i in -14..14 {
            let offset = i as Real / 10.0;
            // Some slices will fail to hit the model.
            if let Ok(slice) = compute_slice(&model, offset) {
                let slice_mesh = geo_bevy::multi_polygon_to_mesh(&slice).unwrap();

                for mesh in slice_mesh {
                    let mut mesh = mesh.mesh;
                    if let VertexAttributeValues::Float32x3(values) =
                        mesh.attribute_mut(Mesh::ATTRIBUTE_NORMAL).unwrap()
                    {
                        values.iter_mut().for_each(|value| *value = [0.0, 0.0, 1.0]);
                    }

                    cmds.spawn((
                        Mesh3d(meshes.add(mesh).clone()),
                        MeshMaterial3d(mats.add(StandardMaterial { ..default() })),
                        Transform::from_translation(Vec3::new(0.0, 0.0, 3.0 + offset as f32)),
                        Wireframe,
                        WireframeColor {
                            color: BLACK.into(),
                        },
                    ));
                }
            }
        }
    }

    cmds.spawn((PointLight::default(), Transform::from_xyz(2., 5., 2.)));
    cmds.spawn((
        Transform::from_translation(Vec3::new(0., 2., 3.)),
        PanOrbitCamera::default(),
    ));
}
