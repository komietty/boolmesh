# boolmesh

[![License: MPL-2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Crates.io](https://img.shields.io/crates/v/boolmesh.svg)](https://crates.io/crates/boolmesh)
[![Docs.rs](https://img.shields.io/docsrs/boolmesh)](https://docs.rs/boolmesh)
[![Test](https://github.com/komietty/boolmesh/actions/workflows/test.yml/badge.svg)](https://github.com/komietty/boolmesh/actions/workflows/test.yml)

![ms](https://raw.githubusercontent.com/komietty/boolmesh/v0.2.0/examples/menger_sponge.png)
Boolmesh is a pure Rust library for performing robust and efficient mesh boolean operations. It is a full-from-scratch Rust implementation inspired by [Elalish’s Manifold](https://manifoldcad.org/docs/html/classmanifold_1_1_manifold.html), which is well known for its robustness and is now part of OpenSCAD.

The codebase is clean and minimal, with dependencies only on `glam`—which provides SIMD acceleration—and optionally on `rayon` for multi-threading support. Besides being robust, Boolmesh is also very fast. For example, generating a Menger Sponge of depth 4 (the model shown above) takes only around 8 seconds on an Apple Silicon M4 using single-threading, and around 4 seconds with multi-threading enabled. Also, the library supports both 32-bit and 64-bit floating point arithmetic.

## Carry any data through booleans

The useful point of boolmesh is that you can attach any value to a mesh and it survives the boolean. UVs, colors, normals, matrices, or your own struct: anything that can be added and scaled by a real number (`Clone + PartialEq + Default + Add + Mul<Real>`) works, since `Manifold<S>` is generic over the value type. Vertices created along the cut get values interpolated from the original faces, and seams such as UV borders stay sharp.

![uv](https://raw.githubusercontent.com/komietty/boolmesh/v0.2.0/examples/value_inheritance.gif)


## Usage
The usage is intentionally simple, as the library exposes only one main function for end users. To perform a boolean operation, construct a mesh buffer structure (called a `Manifold`) from vertex positions and face indices, then call `compute_boolean()` to obtain the result.

Note: Input meshes must be manifold, meaning they must not contain boundaries or overlapping geometry.

``` rust  
// plain meshes
let mfd_0 = Manifold::new(&positions_0, &indices_0).unwrap();
let mfd_1 = Manifold::new(&positions_1, &indices_1).unwrap();

// or meshes carrying one value per triangle corner (uv, color, normal, ...): values_0[i] belongs to indices_0[i]
let mfd_0 = Manifold::new_with_values(&positions_0, &indices_0, &values_0).unwrap();
let mfd_1 = Manifold::new_with_values(&positions_1, &indices_1, &values_1).unwrap();

let result = compute_boolean(&mfd_0, &mfd_1, OpType::Subtract).unwrap();

// extract positions, triangle indices and values
let positions = &result.ps;
let indices = result.hs.iter().map(|h| h.tail).collect::<Vec<_>>();
let values = &result.vals; // one value per triangle corner: values[i] belongs to indices[i]
```

Example of a Menger Sponge generator can be found in the examples folder.

```
 cargo run --package boolmesh --release --example menger_sponge --features=bevy,rayon,f32
```