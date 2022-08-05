#include <vector>
#include <tuple>
#include <numeric>
#include <cmath>
#include <functional>
#include <Eigen/Dense>
#include <boost/geometry.hpp>
#include <boost/range/algorithm/nth_element.hpp>
#include <iostream>
#include <limits>

constexpr double PI = 3.141592653589793;
constexpr double TAU = 6.283185307179586;

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
    P_WithinCone,
    P_WithinConeFrustum,
    P_WithinBall,
    // P_Nearest
};

enum Aggregation
{
    A_Nearest,
    // A_Any,
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

    /// The origin of the RTree
    point_t _origin;
public:
    SphericalRTree() {} // default constructor for Cython wrapper

    /// Create an instance from a vector of points, the values will be filled as the index of the points
    SphericalRTree(std::vector<point_t> points) :
        _points(Eigen::Map<Eigen::Matrix<float, -1, 3, Eigen::RowMajor>>((float*)points.data(), points.size(), 3)),
        _values(points.size()) {
        std::iota(_values.begin(), _values.end(), 0);
    }
    /// Create an instance from a vector of points with corresponding values, the values must have the same size
    /// and the points.
    SphericalRTree(std::vector<point_t> points, std::vector<data_t> values) :
        _points(Eigen::Map<Eigen::Matrix<float, -1, 3, Eigen::RowMajor>>((float*)points.data(), points.size(), 3)),
        _values(std::move(values)) {}

    /// Intialize the RTree with an origin position, all the points will be converted
    /// to a spherical coordinate centered on this origin.
    void initialize(point_t origin, float max_distance = 0.0) {
        _origin = origin;

        // convert points to spherical coordinate
        auto dx = _points.col(0).array() - bg::get<0>(origin);
        auto dy = _points.col(1).array() - bg::get<1>(origin);
        auto dz = _points.col(2).array() - bg::get<2>(origin);

        auto atan2 = [](float y, float x) { return std::atan2(y,x); };
        auto phi = dy.binaryExpr(dx, atan2);
        Eigen::ArrayXf phi_2pi = (phi >= 0).select(phi, phi + TAU);
        Eigen::ArrayXf rxy2 = dx.square() + dy.square();
        Eigen::ArrayXf theta = dz.binaryExpr(rxy2.sqrt(), atan2);
        Eigen::ArrayXf r = (rxy2 + dz.square()).sqrt();

        // construct the RTree with the spherical coordinates
        std::vector<srtree_value_t> projections;
        projections.reserve(_points.rows());
        
        for (size_t i = 0; i < _points.rows(); i++) {
            if (max_distance > 0 && r(i) > max_distance) {
                continue;
            }
            spoint_t sp = { phi_2pi(i), theta(i), r(i) };
            projections.push_back(std::make_pair(sp, _values[i]));
        }
        _tree = srtree_t(projections.begin(), projections.end());
    }

    /// Dispose the existing spherical RTree to free the occupied memory
    void dispose() {
        _tree.clear();
    }
    /// Check whether the RTree is initialized
    bool is_initialized() const {
        return _tree.size() > 0;
    }

    size_t size() const {
        return _points.rows();
    }
    point_t point_at(size_t index) const {
        return point_t { _points(index, 0), _points(index, 1), _points(index, 2) };
    }
    unsigned int value_at(size_t index) const {
        return _values[index];
    }
    point_t origin() const {
        return _origin;
    }

    // position -> values
    std::vector<data_t> query(point_t point, Predicates predicate, float predicate_arg) const {
        float x = bg::get<0>(point) - bg::get<0>(_origin);
        float y = bg::get<1>(point) - bg::get<1>(_origin);
        float z = bg::get<2>(point) - bg::get<2>(_origin);
        std::vector<data_t> indices;

        // cartesian to float
        float phi = std::atan2(y, x);
        phi = phi >= 0 ? phi : phi + TAU;
        float rxy2 = x * x + y * y;
        float theta = std::atan2(z, std::sqrt(rxy2));
        float r = std::sqrt(rxy2 + z * z);

        if (predicate == P_WithinCone || predicate == P_WithinConeFrustum || predicate == P_WithinBall) {
            float range = predicate_arg;

            float rs = range / r;
            float rs_phi = rs / cos(theta);
            float rmin, rmax;

            switch (predicate) {
                case P_WithinCone:
                    rmin = 0; rmax = r + range;
                    break;
                case P_WithinConeFrustum:
                    rmin = r - range; rmax = std::numeric_limits<float>::infinity();
                    break;
                case P_WithinBall:
                    rmin = r - range; rmax = r + range;
                    break;
            }

            spoint_t p1 = {phi - rs_phi, std::max(theta - rs, -float(PI / 2)), rmin};
            spoint_t p2 = {phi + rs_phi, std::min(theta + rs, float(PI / 2)), rmax};
            sbox_t bq(p1, p2);

            auto rit = _tree.qbegin(bgi::within(bq));
            for (auto rit = _tree.qbegin(bgi::within(bq)); rit != _tree.qend(); ++rit)
                indices.push_back(rit->second);
        }
        return indices;
    }

    // list of positions -> values
    std::vector<data_t> query_batch(std::vector<point_t> points,
        Predicates predicate, std::vector<float> predicate_args,
        Aggregation aggregation) const;
};