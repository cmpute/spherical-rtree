
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED NO_MODULE)

add_cython_target(_binding CXX)
add_library(_binding MODULE ${_binding})
target_include_directories(_binding PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
target_link_libraries(_binding ${Boost_LIBRARIES})
python_extension_module(_binding)
install(TARGETS _binding LIBRARY DESTINATION sphertree)
