from libcpp cimport bool
from libcpp.vector cimport vector
from numpy cimport ndarray
from warnings import warn

cdef extern from "SphericalRTree.hpp":
    cdef cppclass point_t:
        pass

    cdef enum Predicates:
        P_WithinCone
        P_WithinConeFrustum
        P_WithinBall

    point_t make_point(float, float, float)
    void unwrap_point(const point_t&, float&, float&, float&)

    cdef cppclass SphericalRTree:
        SphericalRTree()
        SphericalRTree(vector[point_t])
        SphericalRTree(vector[point_t], vector[unsigned int])
        
        size_t size()
        point_t point_at(size_t)
        unsigned int value_at(size_t)
        point_t origin()

        void initialize(point_t origin, float max_distance)
        void dispose()
        bool is_initialized()

        vector[unsigned int] query(point_t point, Predicates predicate, float predicate_args)

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

    def __len__(self):
        return self._tree.size()

    def __getitem__(self, size_t index):
        cdef float x = 0, y = 0, z = 0
        unwrap_point(self._tree.point_at(index), x, y, z)
        return (x, y, z, self._tree.value_at(index))

    def initialize(self, origin, float max_distance = 0.0):
        x, y, z = origin
        cdef point_t c_origin = make_point(x, y, z)
        self._tree.initialize(c_origin, max_distance)

    def dispose(self):
        self._tree.dispose()

    def is_initialized(self):
        self._tree.is_initialized()

    def origin(self):
        cdef float x = 0, y = 0, z = 0
        if self._tree.is_initialized():
            unwrap_point(self._tree.origin(), x, y, z)
            return (x, y, z)
        else:
            return None

    def query(self, point, pred, pred_arg):
        x, y, z = point

        if x == 0 and y == 0:
            warn("the case where x = y = 0 is not currently supported.")

        cdef point_t c_pt = make_point(x, y, z)
        cdef Predicates c_pred = P_WithinCone
        if pred == "cone":
            c_pred = P_WithinCone
        elif pred == "cone_frustum":
            c_pred = P_WithinConeFrustum
        elif pred == "ball":
            c_pred = P_WithinBall
        else:
            raise ValueError("Not supported")

        cdef vector[unsigned int] result = self._tree.query(c_pt, c_pred, pred_arg)
        return result
