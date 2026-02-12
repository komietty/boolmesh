//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use bevy::prelude::*;
use bevy::pbr::wireframe::{WireframePlugin, Wireframe, WireframeColor};
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::PrimitiveTopology;
use bevy::color::palettes::css::*;
use bevy_panorbit_camera::{PanOrbitCamera, PanOrbitCameraPlugin};
use boolmesh::AttrManifold;
use boolmesh::prelude::*;

#[derive(Component)]
struct ToggleableMesh;

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
    mut meshes: ResMut<Assets<Mesh>>
) {
    let mut mfd0 = generate_cube().unwrap();
    let mut mfd1 = generate_cylinder(0.4, 1., 50, 10).unwrap();
    let mut mfd2 = generate_cube().unwrap();
    let idcs1 = mfd1.hs.chunks(3).map(|cs| [cs[0].tail, cs[1].tail, cs[2].tail]).collect::<Vec<_>>();

    let attr0 = (0..mfd0.nf).map(|i| Vec2::new(1., i as f32)).collect::<Vec<_>>();
    let attr1 = (0..mfd1.nf).map(|i| {
        let idcs = idcs1[i];
        let p0 = mfd1.ps[idcs[0]];
        let p1 = mfd1.ps[idcs[1]];
        let p2 = mfd1.ps[idcs[2]];
        Vec2::new(2., (p0 + p1 + p2).y as f32 / 3. + 0.5)
    }).collect::<Vec<_>>();
    let attr2 = (0..mfd2.nf).map(|i| Vec2::new(3., i as f32)).collect::<Vec<_>>();


    mfd1.translate(0.5, 0., 0.25);
    mfd2.translate(-0.5, 0.5, 0.5);

    let amfd0 = AttrManifold::new(mfd0, attr0, true);
    let amfd1 = AttrManifold::new(mfd1, attr1, true);
    let amfd2 = AttrManifold::new(mfd2, attr2, true);

    let res = compute_boolean_with_attributes(&amfd0, &amfd1, OpType::Subtract).unwrap();
    //let res = compute_boolean_with_attributes(&res  , &amfd2, OpType::Subtract).unwrap();

    cmds.spawn((DirectionalLight::default(), Transform::from_xyz(30., 40., 30.)));
    cmds.spawn((Transform::from_translation(Vec3::new(0., 0., 2.)), PanOrbitCamera::default(),));

    let mut m = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::default());

    let mut pos = vec![];
    let mut col = vec![];
    let mut vns = vec![];
    for (fid, hs) in res.manifold.hs.chunks(3).enumerate() {
        let p0  = res.manifold.ps[hs[0].tail];
        let p1  = res.manifold.ps[hs[1].tail];
        let p2  = res.manifold.ps[hs[2].tail];
        let n   = res.manifold.face_normals[fid];
        let var = res.attribute[fid];
        let mut c = [0., 0., 0., 1.];
        if var.x as usize == 1 { c[0] = 1.; }
        if var.x as usize == 2 { c[1] = var.y as f32; }
        if var.x as usize == 3 { c[2] = 1.; }

        pos.push([p0.x as f32, p0.y as f32, p0.z as f32]);
        pos.push([p1.x as f32, p1.y as f32, p1.z as f32]);
        pos.push([p2.x as f32, p2.y as f32, p2.z as f32]);
        col.push(c);
        col.push(c);
        col.push(c);
        vns.push([n.x as f32, n.y as f32, n.z as f32]);
        vns.push([n.x as f32, n.y as f32, n.z as f32]);
        vns.push([n.x as f32, n.y as f32, n.z as f32]);
    }
    m.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    m.insert_attribute(Mesh::ATTRIBUTE_NORMAL, vns);
    m.insert_attribute(Mesh::ATTRIBUTE_COLOR, col);

    cmds.spawn((
        Mesh3d(meshes.add(m).clone()),
        MeshMaterial3d(mats.add(StandardMaterial { unlit: true, base_color: GRAY.into(), ..default() })),
        Transform::default(),
        Wireframe,
        WireframeColor { color: Srgba::rgb(0.3, 0.3, 0.3).into() },
        ToggleableMesh,
    ));

}

