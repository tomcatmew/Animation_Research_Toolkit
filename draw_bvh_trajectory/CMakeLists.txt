  
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

project(draw_bvh_trajectory)

#glfw
IF(MSVC)
  set(GLFW_LIBRARY    ../glfw/src/Release/glfw3.lib)
ELSE()
  set(GLFW_LIBRARY    ../glfw/src/libglfw3.a)
ENDIF()
get_filename_component(GLFW_LIBRARY ${GLFW_LIBRARY} ABSOLUTE)
get_filename_component(GLFW_INCLUDE_DIR ../glfw/include ABSOLUTE)

#opengl
find_package(OpenGL REQUIRED)

#myutility 
set(UTI_INCLUDE_DIR "../utility")
set(SOURCES ${UTI_INCLUDE_DIR}/ycmatrix.cpp ${UTI_INCLUDE_DIR}/ycdraw.cpp ${UTI_INCLUDE_DIR}/ycmatrix.h ${UTI_INCLUDE_DIR}/ycdraw.h)

######################################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${UTI_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
)
add_executable(untilities ${SOURCES})
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