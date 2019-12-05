# libigl + RectangleBinPack

This project is base on the [libigl-example-project](https://github.com/libigl/libigl-example-project) and uses a sligtly modified(CMake integration) version of https://github.com/juj/RectangleBinPack as a submodule.

## What it does

1) read the 3D model
2) cut it along the mean X,Y and Z coordinates to have a vector of patches
3) run LSCM per patch
4) pack these parameterized patches using [RectangleBinPack](https://github.com/BruegelN/RectangleBinPack)


## See the tutorial first

Then build, run and understand the [libigl
tutorial](http://libigl.github.io/libigl/tutorial/) and [example project](https://libigl.github.io/example-project/).

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

[RectangleBinPack with CMake](https://github.com/BruegelN/RectangleBinPack) is included as a submodule

## Get Started

```bash
git clone --recurse-submodules https://github.com/BruegelN/libigl-RectangleBinPack
```
```bash
cd libigl-RectangleBinPack
mkdir build
cd build
cmake ..
make
```

This should find and build the dependencies and create a `example_bin` binary.
You need a local copy of [libigl](https://github.com/libigl/libigl).
For futher information on how add libigl see [the documentation of the libigl example project](https://libigl.github.io/example-project/).

# Run

When execute `./example_bin ../lilium.obj` you can view the n-patch [1-9] by pressing the corresponding number key on your keyboard. Pressing __0__ will show all packed patches. With ' ' (space) you can view the original mesh.
![Packing after cuttuing along X,Y,Z mean values and LSCM Parameterization](https://raw.githubusercontent.com/BruegelN/libigl-RectangleBinPack/master/.github/screenshot.png)
