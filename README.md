# libs
Dependencies folder

This repository contains dependencies for other C++ projects.

## Internal dependencies

* QGL_toolkit: a simplified version of [LibQGLViewer](https://github.com/GillesDebunne/libQGLViewer) with a reduced set of essential features, re-implemented as a small header-only library.
  
  Requires Qt 6.7.0, using the [Qt open source](https://www.qt.io/download-open-source) installer for Windows

* GLtools: minimalist toolkit containing classes for camera, trackball, logger ...*etc*, in a single header file

## External dependencies

All external dependencies are open-source libraries. Latest versions can be installed using the links and intructions listed below. 
For convenience, they are already included in the *third_party* directory. These libraries are already included in the Cmake files provided with each project.

* [GLEW (The OpenGL Extension Wrangler Library)](http://glew.sourceforge.net/)

  Download binaries for windows
  
* [GLFW (Graphics Library Framework)](https://www.glfw.org/)

  Download latest version. The library is provided with a Cmake file, use it to generate a VisualStudio solution, then open it, and build
  
* [GLM (OpenGL Mathematics)](https://github.com/g-truc/glm)

  Just download the libary (header only)

* [Dear ImGui (Immediate-mode Graphical User Interface)](https://github.com/ocornut/imgui)

  Include sources in your project, with the appropriate backend for your configuration
  
* [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/)

  Winged-edge data structure for mesh processing.
  Built from [OpenMesh 11.0.0 sources](https://www.graphics.rwth-aachen.de/software/openmesh/download/). Use the provided Cmake script to generate a VS2022 solution (uncheck the *BUILD_APPS* option), and build the *OpenMeshCore* and *OpenMeshTools* projects to generate the corresponding lib files. 

* [Eigen](https://gitlab.com/libeigen/eigen)

  Linear algebra header-only library, includes matrices, vectors, numerical solvers, and related algorithms

* [stb library](https://github.com/nothings/stb) for image loading:

  Just download and include the stb_image.h header
  
* [tinyobjloader](https://github.com/syoyo/tinyobjloader)

  Just download and include the tiny_obj_loader.h header

* [LodePNG](https://lodev.org/lodepng/)  

  https://github.com/lvandeve/lodepng


