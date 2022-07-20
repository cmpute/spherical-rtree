from libcpp.vector cimport vector
from numpy cimport ndarray

cdef extern from "SphericalRTree.hpp":
    cdef cppclass point_t:
        pass

    point_t make_point(float, float, float)
    void unwrap_point(const point_t&, float&, float&, float&)

    cdef cppclass SphericalRTree:
        SphericalRTree()
        SphericalRTree(vector[point_t])
        SphericalRTree(vector[point_t], vector[unsigned int])
        
        point_t point_at(size_t)
        unsigned int value_at(size_t)

cdef class Tree:
    cdef SphericalRTree _tree

    def __init__(self, float[:,:] points, unsigned int[:] values = None):
        cdef vector[point_t] vec_points
        cdef vector[unsigned int] vec_values

        # convert points to cpp vector
        vec_points.reserve(len(points))
        for i in range(len(points)):
            vec_points.push_back(make_point(points[i, 0], points[i, 1], points[i, 2]))

        if values is None:
            # initialize with points only
            self._tree = SphericalRTree(vec_points)
        elif len(values) != len(points):
            raise ValueError("The length of values should be equal to the length of points")
        else:
            # convert values to cpp vector
            vec_values.reserve(len(values))
            for i in range(len(values)):
                vec_values.push_back(values[i])

            # initialize with points and values
            self._tree = SphericalRTree(vec_points, vec_values)

    def __getitem__(self, size_t index):
        cdef float x = 0, y = 0, z = 0
        unwrap_point(self._tree.point_at(index), x, y, z)
        return (x, y, z, self._tree.value_at(index))
