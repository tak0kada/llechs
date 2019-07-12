#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/detail/tolerance_manip.hpp>
#include <iostream>

#include "unitsp.hpp"
#include "cnthd/mesh.hpp"
#include "cnthd/structure.hpp"

using namespace cnthd;

BOOST_AUTO_TEST_CASE(test_find_center)
{
    {
        // Mesh mesh = read_obj("./cnthd/test/data/cube_small.obj");
        Mesh mesh = read_obj("./test/data/cube_small.obj");
        // now center of the mesh is (0, 0, 0)
        mesh.move(Vector{-0.5, -0.5, -0.5});

        Vector center{find_center(mesh)};
        std::cout << center << std::endl;

        BOOST_TEST(center.x == 0., boost::test_tools::tolerance(1e-8));
        BOOST_TEST(center.y == 0., boost::test_tools::tolerance(1e-8));
        BOOST_TEST(center.z == 0., boost::test_tools::tolerance(1e-8));
    }

    {
        // Mesh mesh = read_obj("./cnthd/test/data/sphere.obj");
        Mesh mesh = read_obj("./test/data/sphere.obj");

        Vector center{find_center(mesh)};
        std::cout << center << std::endl;

        BOOST_TEST(center.x == 0., boost::test_tools::tolerance(1e-8));
        BOOST_TEST(center.y == 0., boost::test_tools::tolerance(1e-8));
        BOOST_TEST(center.z == 0., boost::test_tools::tolerance(1e-8));
    }

    {
        // Mesh mesh = read_obj("./cnthd/test/data/sphere.obj");
        Mesh mesh = read_obj("./test/data/DDGSpring2016/sphere_large.obj");

        Vector center{find_center(mesh)};
        std::cout << center << std::endl;

        BOOST_TEST(center.x == 0., boost::test_tools::tolerance(1e-8));
        BOOST_TEST(center.y == 0., boost::test_tools::tolerance(1e-8));
        BOOST_TEST(center.z == 0., boost::test_tools::tolerance(1e-8));
    }

    {
        // Mesh mesh = read_obj("./cnthd/test/data/sphere.obj");
        Mesh mesh = read_obj("./test/data/cell00038.obj20.obj");

        Vector center{find_center(mesh)};
        std::cout << "circumcenter: " << center << std::endl;

        // the center of mass
        Vector g{0,0,0};
        const std::size_t N{mesh.vertex.size()};
        for (const Vertex& v: mesh.vertex)
        {
            g.x += v.x() / N;
            g.y += v.y() / N;
            g.z += v.z() / N;
        }
        std::cout << "center of mass: " << g << std::endl;

    }
}
