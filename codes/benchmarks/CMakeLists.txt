cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(slam-book-benchmarks)

set(CMAKE_CXX_STANDARD 11)
set(PROJECT_VERSION 1.0.0)

find_package(OpenGL REQUIRED)
if(OpenGL_FOUND)
	message("OpenGL found")
	message ("OPENGL_LIBRARIES: " ${OPENGL_LIBRARIES})
endif()

find_package(GLUT REQUIRED)
if(GLUT_FOUND)
	message("GLUT found")
	message ("GLUT_glut_LIBRARY: " ${GLUT_glut_LIBRARY})
endif()

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
	message("EIGEN3 found")
	message ("EIGEN3_INCLUDE_DIR: " ${EIGEN3_INCLUDE_DIR})
	message ("EIGEN3_INCLUDE_DIRS: " ${EIGEN3_INCLUDE_DIRS})
	message ("Eigen3_DIR: " ${Eigen3_DIR})
endif()

include_directories(
	../common/include
	${EIGEN3_INCLUDE_DIR} 
	${EIGEN3_INCLUDE_DIRS} 
	${Eigen3_DIR}
)

add_executable(pose3d_graph_slam
	src/pose3d_graph_slam.cpp
)

target_link_libraries(pose3d_graph_slam
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)







