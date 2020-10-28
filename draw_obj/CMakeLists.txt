  
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

project(assignment1_draw_obj)

configure_file(fragment.fragmentshader fragment.fragmentshader COPYONLY)
configure_file(vertex.vertexshader vertex.vertexshader COPYONLY)
configure_file(fragment.fragmentshader Debug/fragment.fragmentshader COPYONLY)
configure_file(vertex.vertexshader Debug/vertex.vertexshader COPYONLY)
configure_file(bunny.obj bunny.obj COPYONLY)
configure_file(bunny.obj Debug/bunny.obj COPYONLY)

include_directories(${CMAKE_BINARY_DIR}/configuration)

# glfw
IF(MSVC)
  set(GLFW_LIBRARY    ../glfw/src/Release/glfw3.lib)
ELSE()
  set(GLFW_LIBRARY    ../glfw/src/libglfw3.a)
ENDIF()
get_filename_component(GLFW_LIBRARY ${GLFW_LIBRARY} ABSOLUTE)
get_filename_component(GLFW_INCLUDE_DIR ../glfw/include ABSOLUTE)

#assimp
#set(ASSIMP_LIBRARY   ../assimp-3.3.1/build/code/Debug/assimp-vc140-mt.lib )
#get_filename_component(ASSIMP_LIBRARY ${ASSIMP_LIBRARY} ABSOLUTE)
#set(ASSIMP_INCLUDE_DIR "../assimp-3.3.1/include")
#get_filename_component(ASSIMP_INCLUDE_DIR ${ASSIMP_INCLUDE_DIR} ABSOLUTE)

# headers
#set(HEADER_INCLUDE_DIR headers)
#get_filename_component(HEADER_INCLUDE_DIR ${HEADER_INCLUDE_DIR} ABSOLUTE)

#loader
set(LOADER_INCLUDE_DIR loader)
get_filename_component(LOADER_INCLUDE_DIR ${LOADER_INCLUDE_DIR} ABSOLUTE)

#glad
set(GLAD_INCLUDE_DIR "../glad3/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

#glm
set(GLM_INCLUDE_DIR "../glm")
get_filename_component(GLM_INCLUDE_DIR ${GLM_INCLUDE_DIR} ABSOLUTE)

#glew
#set(GLEW_LIBRARY   ../glew-1.13.0/build/Debug/GLEW_1130.lib)
#get_filename_component(GLEW_LIBRARY ${GLEW_LIBRARY} ABSOLUTE)
#set(GLEW_INCLUDE_DIR "../glew-1.13.0/include")
#get_filename_component(GLEW_INCLUDE_DIR ${GLEW_INCLUDE_DIR} ABSOLUTE)

# required packages opengl/assimp
find_package(OpenGL REQUIRED)
message(STATUS "Found GLFW3 in ${GLFW3_INCLUDE_DIR}")

add_definitions(
	-DGLEW_STATIC
)

######################################################
include_directories (${INCLUDE_DIR}) 
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
  loader/shader.cpp
  loader/shader.hpp
  loader/objloader.cpp
  loader/objloader.hpp
  loader/controls.cpp
  loader/controls.hpp
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