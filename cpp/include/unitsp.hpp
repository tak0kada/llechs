#pragma once

#include <cmath>
#include <exception>
#include "cnthd/mesh.hpp"
#include "cnthd/structure.hpp"
#include "cnthd/util.hpp"

namespace cnthd
{
    // https://www.ipsj.or.jp/07editj/promenade/4309.pdf
    Vector find_center(const Mesh& mesh)
    {
        auto norm = [](const Vector& pos, const Vertex& vert) -> real_t
        {
            const real_t dx{vert.x() - pos.x};
            const real_t dy{vert.y() - pos.y};
            const real_t dz{vert.z() - pos.z};
            return dx * dx + dy * dy + dz * dz;
        };

        const Vertex start{mesh.vertex[0]};
        // vector to represent currently estimated sphere center
        Vector p{start.x(), start.y(), start.z()};

        real_t move{0.5};
        while (move > 1e-8)
        {
            for (int t = 0; t < 50; ++t)
            {
                // index to farthest vertex from p
                std::size_t k{0};
                real_t max{0};
                for (std::size_t i = 0; i < mesh.vertex.size(); ++i)
                {
                    if (norm(p, mesh.vertex[i]) > max)
                    {
                        max = norm(p, mesh.vertex[i]);
                        k = i;
                    }
                }
                p.x += (mesh.vertex[k].x() - p.x) * move;
                p.y += (mesh.vertex[k].y() - p.y) * move;
                p.z += (mesh.vertex[k].z() - p.z) * move;
            }
            move *= 0.5;
        }

        return p;
    }

    void to_unit_sphere(Mesh& mesh)
    {
        const auto center = find_center(mesh);
        mesh.move(-center);

        for (Vertex& v: mesh.vertex)
        {
            cnthd::real_t l = std::sqrt(v.p[0] * v.p[0] + v.p[1] * v.p[1] + v.p[2] * v.p[2]);
            if (l == 0)
            {
                throw std::runtime_error("ERROR: There is a vector on the origin.");
            }
            else
            {
                v.p[0] /= l;
                v.p[1] /= l;
                v.p[2] /= l;
            }
        }
    }

    // real_t sph_area(const Face& f)
    // {
    //     const auto [vert0, vert1, vert2] = f.vertex();
    //     const Vector v0{vert0->x(), vert0->y(), vert0->z()};
    //     const Vector v1{vert1->x(), vert1->y(), vert1->z()};
    //     const Vector v2{vert2->x(), vert2->y(), vert2->z()};
    //
    //     // unit normal to each Vector vn
    //     const Vector n02{unit_vec(v2 - v0 * (v0*v2 / v0.norm()))};
    //     const Vector n01{unit_vec(v1 - v0 * (v0*v1 / v0.norm()))};
    //     const Vector n10{unit_vec(v0 - v1 * (v1*v0 / v1.norm()))};
    //     const Vector n12{unit_vec(v2 - v1 * (v1*v2 / v1.norm()))};
    //     const Vector n21{unit_vec(v1 - v2 * (v2*v1 / v2.norm()))};
    //     const Vector n20{unit_vec(v0 - v2 * (v2*v0 / v2.norm()))};
    //
    //     return std::acos(n02 * n01) + std::acos(n10 * n12) + std::acos(n21 * n20) - pi;
    // }
    //
    // real_t theta(const Vertex& v)
    // {
    //     const real_t x{v.x()};
    //     const real_t y{v.y()};
    //     const real_t z{v.z()};
    //     const real_t l{std::sqrt(x*x + y*y + z*z)};
    //     if (l == 0)
    //     {
    //         return 0;
    //     }
    //     else
    //     {
    //         return std::acos(z / l);
    //     }
    // }
    //
    // real_t phi(const Vertex& v)
    // {
    //     const real_t x{v.x()};
    //     const real_t y{v.y()};
    //     if (x == 0 && y == 0)
    //     {
    //         return 0;
    //     }
    //     else if (y > 0)
    //     {
    //         return std::atan2(y, x);
    //     }
    //     else if (y < 0)
    //     {
    //         return 2 * pi + std::atan2(y, x);
    //     }
    //     else
    //     {
    //         return pi;
    //     }
    // }

}
