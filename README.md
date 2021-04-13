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

build project bvh_trajectory_nn folder by
```
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/absolute/path/to/libtorch ..

and copy all .dll files from libtorch to project directory
```


## To Do :
- [ ] understand the paper fourier feature network
- [ ] checkout related paper from phased functioned neural network to find out previous simple model of training bvh files
- [ ] Successfully use the trained netwokr in C++ environment and display a decent animation result.
- [x] Further improve the NN to allow prediction of rotation besides transformation. may need quaternion ?
- [ ] Finish the the training

## Update Note:
2021.4.11. update the test_nn for testing the C++ pytroch inplementation
2021.4.13. The NN successfully predict the XYZ transformation of the root position. Now works on rotation.
