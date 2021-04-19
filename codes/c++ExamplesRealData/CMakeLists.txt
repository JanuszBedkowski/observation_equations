cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(slam-book-examples)

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

find_package( Eigen3 REQUIRED)
if(EIGEN3_FOUND)
	message("EIGEN3 found")
	message ("EIGEN3_INCLUDE_DIR: " ${EIGEN3_INCLUDE_DIR})
	message ("EIGEN3_INCLUDE_DIRS: " ${EIGEN3_INCLUDE_DIRS})
	message ("Eigen3_DIR: " ${Eigen3_DIR})
endif()

find_package(PCL 1.5 REQUIRED)


set(Boost_USE_STATIC_LIBS        OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost REQUIRED COMPONENTS system filesystem thread iostreams) 
if(Boost_FOUND)
    message ("Boost_LIBRARIES: " ${Boost_LIBRARIES})
    message ("Boost_INCLUDE_DIRS: " ${Boost_INCLUDE_DIRS})
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

include_directories(
	${EIGEN3_INCLUDE_DIR} 
	${EIGEN3_INCLUDE_DIRS} 
	${Eigen3_DIR}
	${Boost_INCLUDE_DIRS} 
	${PCL_INCLUDE_DIRS}
)

find_package( OpenCV REQUIRED )

include_directories(
	${EIGEN3_INCLUDE_DIR} 
	${EIGEN3_INCLUDE_DIRS} 
	${Eigen3_DIR}
	${OpenCV_INCLUDE_DIRS}
)

add_executable(particle_filter_demo
	src/particle_filter_demo.cpp
)

target_link_libraries(particle_filter_demo
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
	${OpenCV_LIBS}
)

add_executable(photogrammetry
	src/photogrammetry.cpp
)

target_link_libraries(photogrammetry
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
	${OpenCV_LIBS}
)

add_executable(point_cloud_registration
	src/point_cloud_registration.cpp
)

target_link_libraries(point_cloud_registration
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
	${OpenCV_LIBS}
	${PCL_LIBRARIES}
)

