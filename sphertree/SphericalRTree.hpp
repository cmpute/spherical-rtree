#include <vector>
#include <tuple>
#include <numeric>
#include <Eigen/Core>
#include <boost/geometry.hpp>
#include <boost/range/algorithm/nth_element.hpp>
#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 3, bg::cs::cartesian> point_t; // point in cartesian coordinate
typedef bg::cs::spherical_equatorial<bg::radian> scoord_t;
typedef bg::model::point<float, 3, scoord_t> spoint_t; // point in spherical coordinate
typedef unsigned int data_t; // store id as data
typedef std::pair<spoint_t, data_t> srtree_value_t;
typedef bg::model::box<spoint_t> sbox_t;
typedef bgi::rstar<16, 4> srtree_params_t;
typedef bgi::rtree<srtree_value_t, srtree_params_t> srtree_t; // rtree in spherical coordinate

enum Predicates
{
    Within,
};

enum Aggregation
{
    Closest,
    Any,
};

point_t make_point(float x, float y, float z) {
    return point_t {x, y, z};
}

void unwrap_point(const point_t &pt, float &x, float &y, float &z) {
    x = bg::get<0>(pt);
    y = bg::get<1>(pt);
    z = bg::get<2>(pt);
}

class SphericalRTree
{
private:
    /// The spherical rtree data structure
    srtree_t _tree;

    /// The point cloud
    Eigen::Matrix<float, -1, 3> _points;

    /// The associated data on the cloud
    std::vector<data_t> _values;
public:
    SphericalRTree() {} // default constructor for Cython wrapper

    SphericalRTree(std::vector<point_t> points) :
        _points(Eigen::Map<Eigen::Matrix<float, -1, 3, Eigen::RowMajor>>((float*)points.data(), points.size(), 3)),
        _values(points.size()) {
        std::iota(_values.begin(), _values.end(), 0);
    }
    SphericalRTree(std::vector<point_t> points, std::vector<data_t> values) :
        _points(Eigen::Map<Eigen::Matrix<float, -1, 3, Eigen::RowMajor>>((float*)points.data(), points.size(), 3)),
        _values(std::move(values)) {}

    // initialize with point index as data
    // TODO: use iterator to initialize?
    void initialize(point_t origin);
    void dispose();

    point_t point_at(size_t index) const {
        return point_t { _points(index, 0), _points(index, 1), _points(index, 2) };
    }
    unsigned int value_at(size_t index) const {
        return _values[index];
    }

    // position -> values
    std::vector<data_t> query(point_t point, Predicates predicate, std::vector<float> predicate_args) const;

    // list of positions -> values
    std::vector<data_t> query_batch(std::vector<point_t> points,
        Predicates predicate, std::vector<float> predicate_args,
        Aggregation aggregation) const;
};