  
cmake_minimum_required(VERSION 2.8)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

######################################################

project(draw_cuboid)

configure_file(ColorFragmentShader.fragmentshader ColorFragmentShader.fragmentshader COPYONLY)
configure_file(TransformVertexShader.vertexshader TransformVertexShader.vertexshader COPYONLY)

# glfw
IF(MSVC)
  set(GLFW_LIBRARY    ../glfw/src/Release/glfw3.lib)
ELSE()
  set(GLFW_LIBRARY    ../glfw/src/libglfw3.a)
ENDIF()
get_filename_component(GLFW_LIBRARY ${GLFW_LIBRARY} ABSOLUTE)
get_filename_component(GLFW_INCLUDE_DIR ../glfw/include ABSOLUTE)

#loader
set(LOADER_INCLUDE_DIR load_shader)
get_filename_component(LOADER_INCLUDE_DIR ${LOADER_INCLUDE_DIR} ABSOLUTE)

#glm
set(GLM_INCLUDE_DIR "../glm")
get_filename_component(GLM_INCLUDE_DIR ${GLM_INCLUDE_DIR} ABSOLUTE)

#glad
set(GLAD_INCLUDE_DIR "../glad3/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

# opengl
find_package(OpenGL REQUIRED)
message(STATUS "Found GLFW3 in ${GLFW3_INCLUDE_DIR}")
######################################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${LOADER_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${GLM_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ../glad3/include/glad/glad.h 
  ../glad3/src/glad.c
  load_shader/shader.cpp
  load_shader/shader.hpp
  main.cpp
)

#######################################################

IF(APPLE)
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
  find_library( COCOA_LIBRARY Cocoa )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARY}
    ${COCOA_LIBRARY}
    ${IOKIT_LIBRARY}
    ${COREVIDEO_LIBRARY}
  )
ELSEIF(UNIX)
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARY}
    -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
  )
ELSEIF(WIN32)    
  target_link_libraries(${PROJECT_NAME}     
    ${GLFW_LIBRARY}
    ${OPENGL_LIBRARY}
  )
ELSEIF(WIN64)    
  target_link_libraries(${PROJECT_NAME}     
    ${GLFW_LIBRARY}
    ${OPENGL_LIBRARY}
  )
ENDIF()