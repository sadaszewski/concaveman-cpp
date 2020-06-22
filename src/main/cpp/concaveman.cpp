/*
<%
setup_pybind11(cfg)
%>
*/

//
// Author: Stanislaw Adaszewski, 2019
//

#include "concaveman.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#define DEBUG 0

py::array_t<double> ccm(
	py::array_t<double> points_,
	py::array_t<int32_t> hull_points_,
    double concavity, double lengthThreshold) {

	auto points_c = points_.unchecked<2>();
	auto num_points = points_c.shape(0);
    if (points_c.shape(1) != 2) {
        throw std::invalid_argument("points must be 2d");
    }

	auto hull_points_c = hull_points_.unchecked<1>();
	auto num_hull_points = hull_points_c.shape(0);

#if DEBUG
    std::cout << "pyconcaveman2d(), concavity: " << concavity <<
        " lengthThreshold: " << lengthThreshold << std::endl;
#endif

    typedef double T;
    typedef std::array<T, 2> point_type;

    std::vector<point_type> points(num_points);
    for (auto i = 0; i < num_points; i++) {
        points[i] = { points_c(i, 0), points_c(i, 1) };
    }

#if DEBUG
    std::cout << "points:" << std::endl;
    for (auto &p : points)
        std::cout << p[0] << " " << p[1] << std::endl;
#endif

    std::vector<int32_t> hull(num_hull_points);
    for (auto i = 0; i < num_hull_points; i++) {
        hull[i] = hull_points_c(i);
    }

#if DEBUG
    std::cout << "hull:" << std::endl;
    for (auto &i : hull)
        std::cout << i << std::endl;
#endif

    auto concave_points = concaveman<T, 16>(points, hull, concavity, lengthThreshold);

#if DEBUG
    std::cout << "concave_points:" << std::endl;
    for (auto &p : concave_points)
        std::cout << p[0] << " " << p[1] << std::endl;
#endif

    py::array_t<double> result({concave_points.size(), size_t(2)});
    auto result_ptr = result.mutable_unchecked<2>();

    for (size_t i = 0; i < concave_points.size(); i++) {
        result_ptr(i, 0) = concave_points[i][0];
        result_ptr(i, 1) = concave_points[i][1];
    }

    return result;
}

PYBIND11_MODULE(concaveman, m) {
    m.def("concaveman", &ccm);
}
