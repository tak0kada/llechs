#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/detail/tolerance_manip.hpp>
#include <iostream>
#include <algorithm>

#include "util.hpp"
#include "unitsp.hpp"
#include "cnthd/mesh.hpp"
#include "cnthd/structure.hpp"
#include "cnthd/util.hpp"

using namespace cnthd;

BOOST_AUTO_TEST_CASE(test_util)
{
    {
        Mesh mesh = read_obj("./cnthd/test/data/cube_small.obj");
        mesh.move(Vector{-0.5, -0.5, -0.5});
        // map to unit sphere
        to_unit_sphere(mesh);

        for (const auto& f: mesh.face) {
            BOOST_TEST(sph_area(f) == pi / 3, boost::test_tools::tolerance(1e-12));
        }
    }

    {
        Mesh mesh = read_obj("./cnthd/test/data/sphere.obj");

        real_t sum{};
        for (const auto& f: mesh.face) {
            sum += sph_area(f);
        }
        BOOST_TEST(sum == 4 * pi, boost::test_tools::tolerance(1e-12));
    }
}

