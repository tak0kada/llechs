#pragma once

#include <cmath>
#include "cnthd/mesh.hpp"
#include "cnthd/structure.hpp"
#include "cnthd/util.hpp"


namespace cnthd
{

    real_t sph_area(const Face& f)
    {
        const auto [vert0, vert1, vert2] = f.vertex();
        // vector with length 1
        const Vector v0{vert0->x(), vert0->y(), vert0->z()};
        const Vector v1{vert1->x(), vert1->y(), vert1->z()};
        const Vector v2{vert2->x(), vert2->y(), vert2->z()};

        if (v0 == v1 || v1 == v2 || v2 == v0)
        {
            return 0;
        }

        // unit normal to each Vector vn
        // norm of vx == 1
        const Vector n02{unit_vec(v2 - v0 * (v0*v2))};
        const Vector n01{unit_vec(v1 - v0 * (v0*v1))};
        const Vector n10{unit_vec(v0 - v1 * (v1*v0))};
        const Vector n12{unit_vec(v2 - v1 * (v1*v2))};
        const Vector n21{unit_vec(v1 - v2 * (v2*v1))};
        const Vector n20{unit_vec(v0 - v2 * (v2*v0))};

        if(isnan(std::acos(n02 * n01) + std::acos(n10 * n12) + std::acos(n21 * n20) - pi))
        {
            // this happens when the area is too small
            return 0;
        }
        else
        {
            return std::acos(n02 * n01) + std::acos(n10 * n12) + std::acos(n21 * n20) - pi;
        }
    }

    real_t theta(const Vertex& v)
    {
        const real_t x{v.x()};
        const real_t y{v.y()};
        const real_t z{v.z()};
        const real_t l{std::sqrt(x*x + y*y + z*z)};
        if (l == 0)
        {
            return 0;
        }
        else
        {
            return std::acos(z / l);
        }
    }

    real_t phi(const Vertex& v)
    {
        const real_t x{v.x()};
        const real_t y{v.y()};
        if (x == 0 && y == 0)
        {
            return 0;
        }
        else if (y > 0)
        {
            return std::atan2(y, x);
        }
        else if (y < 0)
        {
            return 2 * pi + std::atan2(y, x);
        }
        else
        {
            return pi;
        }
    }

}
