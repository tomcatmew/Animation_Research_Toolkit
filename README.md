### Animation Research Toolkit
This project is based on the assignments given by Prof.Umetatni. 

### Directory 
  * one_triangle - simple C++ project to draw a triangle
  * draw_obj - simple .OBJ parser 
  * draw_cuboid - simple cuboid called by OpenGL libraries
  * draw_arm - simulating animation transformation and rotation by affine matrix 
  * draw_bvh - simple .BVH parser
### Build
build the GLFW before build each of the project code
```
cmake . -A x64
cmake --build . --config Release
```

build each project folder by
```
cmake -S . -B build 
```
