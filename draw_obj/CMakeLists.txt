  
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

# glfw
IF(MSVC)
  set(GLFW_LIBRARY    ../glfw/src/Release/glfw3.lib)
ELSE()
  set(GLFW_LIBRARY    ../glfw/src/libglfw3.a)
ENDIF()
get_filename_component(GLFW_LIBRARY ${GLFW_LIBRARY} ABSOLUTE)
get_filename_component(GLFW_INCLUDE_DIR ../glfw/include ABSOLUTE)

#assimp
set(ASSIMP_LIBRARY   ../assimp-3.3.1/build/code/Debug/assimp-vc140-mt.lib )
get_filename_component(ASSIMP_LIBRARY ${ASSIMP_LIBRARY} ABSOLUTE)
set(ASSIMP_INCLUDE_DIR "../assimp-3.3.1/include")
get_filename_component(ASSIMP_INCLUDE_DIR ${ASSIMP_INCLUDE_DIR} ABSOLUTE)

# headers
set(HEADER_INCLUDE_DIR headers)
get_filename_component(HEADER_INCLUDE_DIR ${HEADER_INCLUDE_DIR} ABSOLUTE)

#glad
set(GLAD_INCLUDE_DIR "../glad4/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

#glm
set(GLM_INCLUDE_DIR "../glm")
get_filename_component(GLM_INCLUDE_DIR ${GLM_INCLUDE_DIR} ABSOLUTE)

# required packages opengl/assimp
find_package(GLM REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")
find_package(OpenGL REQUIRED)
message(STATUS "Found GLFW3 in ${GLFW3_INCLUDE_DIR}")
find_package(ASSIMP REQUIRED)
message(STATUS "Found ASSIMP in ${ASSIMP_INCLUDE_DIR}")
######################################################
include_directories (${INCLUDE_DIR}) 
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${ASSIMP_INCLUDE_DIR}
  ${HEADER_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${GLM_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ../glad4/include/glad/glad.h 
  ../glad4/src/glad.c
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
    ${ASSIMP_LIBRARY}
  )
ELSEIF(UNIX)
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARY}
    ${ASSIMP_LIBRARY}
    -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
  )
ELSEIF(WIN32)    
  target_link_libraries(${PROJECT_NAME}     
    ${GLFW_LIBRARY}
    ${OPENGL_LIBRARY}
    ${ASSIMP_LIBRARY}
  )
ELSEIF(WIN64)    
  target_link_libraries(${PROJECT_NAME}     
    ${GLFW_LIBRARY}
    ${OPENGL_LIBRARY}
    ${ASSIMP_LIBRARY}
  )
ENDIF()