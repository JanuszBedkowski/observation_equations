cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(slam-book-examples)

set(CMAKE_CXX_STANDARD 11)
set(PROJECT_VERSION 1.0.0)

if (WIN32)
	add_compile_definitions(_USE_MATH_DEFINES)
	include_directories(3rdparty/eigen)
else()
	find_package(Eigen3 REQUIRED)
	if(EIGEN3_FOUND)
		message("EIGEN3 found")
		message ("EIGEN3_INCLUDE_DIR: " ${EIGEN3_INCLUDE_DIR})
		message ("EIGEN3_INCLUDE_DIRS: " ${EIGEN3_INCLUDE_DIRS})
		message ("Eigen3_DIR: " ${Eigen3_DIR})

		include_directories(
			${EIGEN3_INCLUDE_DIR} 
			${EIGEN3_INCLUDE_DIRS} 
			${Eigen3_DIR}
		)
	endif()
endif()

find_package(OpenGL REQUIRED)
if(OpenGL_FOUND)
	message("OpenGL found")
	message ("OPENGL_LIBRARIES: " ${OPENGL_LIBRARIES})
endif()

if (WIN32)
	include_directories(3rdparty/freeglut-3.2.2/include/)
	add_subdirectory(3rdparty/freeglut-3.2.2)
else()
	find_package(GLUT REQUIRED)
	if(GLUT_FOUND)
		message("GLUT found")
		message ("GLUT_glut_LIBRARY: " ${GLUT_glut_LIBRARY})
	endif()
endif()

#find_package( OpenCV REQUIRED )

include_directories(
	../common/include
    ../python-scripts/
	../python-scripts/camera-metrics
	../python-scripts/constraints
	../python-scripts/example_func
    ../python-scripts/feature-to-feature-metrics
	../python-scripts/metrics-fusion
	../python-scripts/point-to-feature-metrics
	../python-scripts/point-to-point-metrics
	../python-scripts/slerp
	#${EIGEN3_INCLUDE_DIR} 
	#${EIGEN3_INCLUDE_DIRS} 
	#${Eigen3_DIR}
	#${OpenCV_INCLUDE_DIRS}
)

add_executable(perspective_camera_ba
	src/perspective_camera_ba.cpp
)

if (WIN32)
target_link_libraries(perspective_camera_ba	
	${OPENGL_LIBRARIES} 
	freeglut_static
)
else()
target_link_libraries(perspective_camera_ba
	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(metric_camera_colinearity_ba
	src/metric_camera_colinearity_ba.cpp
)

if (WIN32)
target_link_libraries(metric_camera_colinearity_ba
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(metric_camera_colinearity_ba
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(metric_camera_coplanarity_ba
	src/metric_camera_coplanarity_ba.cpp
)

if (WIN32)
target_link_libraries(metric_camera_coplanarity_ba
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(metric_camera_coplanarity_ba
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(equirectangular_camera_colinearity_ba
	src/equirectangular_camera_colinearity_ba.cpp
)

if (WIN32)
target_link_libraries(equirectangular_camera_colinearity_ba
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(equirectangular_camera_colinearity_ba
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(equirectangular_camera_coplanarity_ba
	src/equirectangular_camera_coplanarity_ba.cpp
)

if (WIN32)
target_link_libraries(equirectangular_camera_coplanarity_ba
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(equirectangular_camera_coplanarity_ba
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_point
	src/point_to_point.cpp
)
if (WIN32)
target_link_libraries(point_to_point
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_point
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_point_with_scale
	src/point_to_point_with_scale.cpp
)
if (WIN32)
target_link_libraries(point_to_point_with_scale
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_point_with_scale
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(perspective_camera_external_orientation_plucker_line
	src/perspective_camera_external_orientation_plucker_line.cpp
)
if (WIN32)
target_link_libraries(perspective_camera_external_orientation_plucker_line
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(perspective_camera_external_orientation_plucker_line
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_point_source_to_target
	src/point_to_point_source_to_target.cpp
)
if (WIN32)
target_link_libraries(point_to_point_source_to_target
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_point_source_to_target
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(plucker_line_to_plucker_line
	src/plucker_line_to_plucker_line.cpp
)
if (WIN32)
target_link_libraries(plucker_line_to_plucker_line
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(plucker_line_to_plucker_line
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(plane_to_plane
	src/plane_to_plane.cpp
)
if (WIN32)
target_link_libraries(plane_to_plane
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(plane_to_plane
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_line
	src/point_to_line.cpp
)
if (WIN32)
target_link_libraries(point_to_line
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_line
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_plane
	src/point_to_plane.cpp
)
if (WIN32)
target_link_libraries(point_to_plane
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_plane
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_projection_onto_line
	src/point_to_projection_onto_line.cpp
)
if (WIN32)
target_link_libraries(point_to_projection_onto_line
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_projection_onto_line
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_projection_onto_plane
	src/point_to_projection_onto_plane.cpp
)
if (WIN32)
target_link_libraries(point_to_projection_onto_plane
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_projection_onto_plane
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(simple_optimization_problem_func_x_with_constraints
	src/simple_optimization_problem_func_x_with_constraints.cpp
)
if (WIN32)
target_link_libraries(simple_optimization_problem_func_x_with_constraints
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(simple_optimization_problem_func_x_with_constraints
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(simple_optimization_problem_func_xy_with_constraints
	src/simple_optimization_problem_func_xy_with_constraints.cpp
)
if (WIN32)
target_link_libraries(simple_optimization_problem_func_xy_with_constraints
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(simple_optimization_problem_func_xy_with_constraints
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(relative_pose
	src/relative_pose.cpp
)
if (WIN32)
target_link_libraries(relative_pose
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(relative_pose
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(relative_pose_covariances
	src/relative_pose_covariances.cpp
)
if (WIN32)
target_link_libraries(relative_pose_covariances
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(relative_pose_covariances
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(smoothness
	src/smoothness.cpp
)
if (WIN32)
target_link_libraries(smoothness
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(smoothness
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(georeference-case1
	src/georeference-case1.cpp
)
if (WIN32)
target_link_libraries(georeference-case1
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(georeference-case1
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(georeference-case2
	src/georeference-case2.cpp
)
if (WIN32)
target_link_libraries(georeference-case2
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(georeference-case2
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(georeference-case4
	src/georeference-case4.cpp
)
if (WIN32)
target_link_libraries(georeference-case4
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(georeference-case4
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(surface_reconstruction_from_lidar_point_cloud
	src/surface_reconstruction_from_lidar_point_cloud.cpp
)
if (WIN32)
target_link_libraries(surface_reconstruction_from_lidar_point_cloud
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(surface_reconstruction_from_lidar_point_cloud
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(optimizing_scale_VO_with_lidar_data
	src/optimizing_scale_VO_with_lidar_data.cpp
)
if (WIN32)
target_link_libraries(optimizing_scale_VO_with_lidar_data
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(optimizing_scale_VO_with_lidar_data
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

#add_executable(perspective_camera_intrinsic_calibration
#	src/perspective_camera_intrinsic_calibration.cpp
#)
#if (WIN32)
#target_link_libraries(perspective_camera_intrinsic_calibration
#  	${OPENGL_LIBRARIES}
#	freeglut_static
#)
#else()
#target_link_libraries(perspective_camera_intrinsic_calibration
#  	${OPENGL_LIBRARIES}
#	${GLUT_glut_LIBRARY}
#)
#endif()

add_executable(bundle_of_rays_intersection
	src/bundle_of_rays_intersection.cpp
)
if (WIN32)
target_link_libraries(bundle_of_rays_intersection
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(bundle_of_rays_intersection
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(rectangular_object_with_unknown_width_height
	src/rectangular_object_with_unknown_width_height.cpp
)
if (WIN32)
target_link_libraries(rectangular_object_with_unknown_width_height
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(rectangular_object_with_unknown_width_height
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(sheaf_of_planes
	src/sheaf_of_planes.cpp
)
if (WIN32)
target_link_libraries(sheaf_of_planes
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(sheaf_of_planes
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(planes_intersection
	src/planes_intersection.cpp
)
if (WIN32)
target_link_libraries(planes_intersection
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(planes_intersection
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(view_nns
	src/view_nns.cpp
)
if (WIN32)
target_link_libraries(view_nns
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(view_nns
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(show_axes
	src/show_axes.cpp
)
if (WIN32)
target_link_libraries(show_axes
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(show_axes
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(distance_point_to_plane
	src/distance_point_to_plane.cpp
)
if (WIN32)
target_link_libraries(distance_point_to_plane
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(distance_point_to_plane
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(m_estimators_functions_viewer
	src/m_estimators_functions_viewer.cpp
)
if (WIN32)
target_link_libraries(m_estimators_functions_viewer
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(m_estimators_functions_viewer
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(adaptive_robust_loss_function_viewer
	src/adaptive_robust_loss_function_viewer.cpp
)
if (WIN32)
target_link_libraries(adaptive_robust_loss_function_viewer
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(adaptive_robust_loss_function_viewer
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(adaptive_robust_loss_function_demo
	src/adaptive_robust_loss_function_demo.cpp
)
if (WIN32)
target_link_libraries(adaptive_robust_loss_function_demo
	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(adaptive_robust_loss_function_demo
	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(point_to_point_source_to_landmark
	src/point_to_point_source_to_landmark.cpp
)
if (WIN32)
target_link_libraries(point_to_point_source_to_landmark
  	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(point_to_point_source_to_landmark
  	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(ceres_robust_loss_function_viewer
	src/ceres_robust_loss_function_viewer.cpp
)
if (WIN32)
target_link_libraries(ceres_robust_loss_function_viewer
   ${OPENGL_LIBRARIES}
   freeglut_static
)
else()
target_link_libraries(ceres_robust_loss_function_viewer
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(distance_to_circle
	src/distance_to_circle.cpp
)
if (WIN32)
target_link_libraries(distance_to_circle
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(distance_to_circle
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(mirror_sim
	src/mirror_sim.cpp
)
if (WIN32)
target_link_libraries(mirror_sim
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(mirror_sim
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

add_executable(pose_interpolation
	src/pose_interpolation.cpp
)
if (WIN32)
target_link_libraries(pose_interpolation
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(pose_interpolation
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()

include_directories(
3rdparty/imgui 
3rdparty/imgui/backends
)
set(IMGUI_SRCS
3rdparty/imgui/imgui.cpp 
3rdparty/imgui/imgui_demo.cpp 
3rdparty/imgui/imgui_draw.cpp 
3rdparty/imgui/imgui_tables.cpp 
3rdparty/imgui/imgui_widgets.cpp
3rdparty/imgui/backends/imgui_impl_glut.cpp 
3rdparty/imgui/backends/imgui_impl_opengl2.cpp
)

add_executable(elementary_error_theory_for_terrestrial_laser_scanner
	src/elementary_error_theory_for_terrestrial_laser_scanner.cpp
	${IMGUI_SRCS}
)
if (WIN32)
target_link_libraries(elementary_error_theory_for_terrestrial_laser_scanner
   	${OPENGL_LIBRARIES}
	freeglut_static
)
else()
target_link_libraries(elementary_error_theory_for_terrestrial_laser_scanner
   	${OPENGL_LIBRARIES}
	${GLUT_glut_LIBRARY}
)
endif()