#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <tuple>
#include <boost/math/constants/constants.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include "cnthd/mesh.hpp"
#include "spherical_harmonics/sph.hpp"
#include "unitsp.hpp"
#include "util.hpp"

namespace cnthd
{


namespace distribution
{
    namespace util
    {
        real_t coef(const Mesh& mesh, const unsigned int& l, const int& m, const std::vector<real_t>& data, const std::vector<real_t>& face_area)
        {
            // the order of vertices in the mesh.vertex and data should be identical
            assert(mesh.nV == data.size());
            // variables for Kahan summation algorithm
            real_t coef{};
            real_t c{};
            volatile real_t y{};
            volatile real_t t{};
            // real_t y{};
            // real_t t{};

            for (const auto& f: mesh.face)
            {
                const std::array<cnthd::Vertex*, 3> vs{
                    std::get<0>(f.vertex()), std::get<1>(f.vertex()), std::get<2>(f.vertex())};

                for (std::size_t i = 0; i < 3; ++i)
                {
                    // Kahan summation algorithm
                    // assume there is a barycentric dual mesh
                    y = sph::sph_harm(l, m, cnthd::theta(*vs[i]), cnthd::phi(*vs[i])) * data[vs[i]->idx] * face_area[f.idx] / 3 - c;
                    t = coef + y;
                    c = (t - coef) - y;
                    coef = t;
                    // end of Kahan algorithm
                }
            }
            return coef;
        }
    }

    std::vector<real_t> calc_coef(const Mesh& mesh, const int& n_coef, const std::vector<real_t>& data)
    {
        // the order of vertices in the mesh.vertex and data should be identical
        assert(mesh.nV == data.size());

        // cache area of faces
        std::vector<real_t> face_area(mesh.nF);
        for (const auto& f: mesh.face)
        {
            face_area[f.idx] = sph_area(f);
        }

        // value to return
        std::vector<real_t> SH(n_coef, 0);

        int l = 0;
        int m = 0;
        for (int n = 0; n < n_coef; ++n)
        {
            SH[n] = util::coef(mesh, l, m, data, face_area);

            if (m == static_cast<int>(l))
            {
                ++l;
                m = -l;
            }
            else
            {
                ++m;
            }
        }

        return SH;
    }
} // end of namespace distribution


namespace shape
{
    struct shmap
    {
    public:
        std::vector<real_t> face_area{};
        std::vector<real_t> x{};
        std::vector<real_t> y{};
        std::vector<real_t> z{};
        Mesh sphere;
        Eigen::Matrix3d Rot{};

        shmap(const Mesh& orig, const Mesh& sphere)
        : sphere{sphere}
        {
            // NOTE: other than the restriction below,
            // the order of vertices in the mesh.vertex and data should be identical
            assert(orig.nE == sphere.nE);
            assert(orig.nF == sphere.nF);
            assert(orig.nV == sphere.nV);

            x.resize(orig.nV);
            y.resize(orig.nV);
            z.resize(orig.nV);
            for (std::size_t i = 0; i < orig.nV; ++i)
            {
                x[i] = orig.vertex[i].x();
                y[i] = orig.vertex[i].y();
                z[i] = orig.vertex[i].z();
            }

            // cache area of faces
            face_area.resize(sphere.nF);
            for (const auto& f: sphere.face)
            {
                face_area[f.idx] = sph_area(f);
            }
        }

        void align_domain();
        std::array<std::vector<real_t>, 3> calc_coef(const unsigned int& n_coef) const;

    private:
        real_t coef(const unsigned int& l, const int& m, const std::vector<real_t>& data) const;
    };

    double shmap::coef(const unsigned int& l, const int& m, const std::vector<real_t>& data) const
    {
        return distribution::util::coef(this->sphere, l, m, data, face_area);
    }

    std::array<std::vector<real_t>, 3> shmap::calc_coef(const unsigned int& n_coef) const
    {
        // value to return
        std::array<std::vector<real_t>, 3> SH{};
        for (std::size_t i = 0; i < 3; ++i)
        {
            SH[i].resize(n_coef);
        }

        const std::array<const std::vector<real_t>* const, 3> xyz{&x, &y, &z};
        for (std::size_t i = 0; i < 3; ++i)
        {
            // use int for l to avoid overflow in calculation of -l
            int l = 0;
            int m = 0;
            for (unsigned int n = 0; n < n_coef; ++n)
            {
                SH[i][n] = coef(l, m, *(xyz[i]));
                if (m == l)
                {
                    ++l;
                    m = -l;
                }
                else
                {
                    ++m;
                }
            }
        }
        return SH;
    }

    void shmap::align_domain()
    {
        // rotate sphere so that we can make alignment by Althloothi's algorithm
        // 1D (l=1) spherical harmonics represent a spheroid
        Eigen::Matrix3d SH1D{};
        SH1D << coef(1, 1, x), coef(1, 1, y), coef(1, 1, z),
                coef(1, -1, x), coef(1, -1, y), coef(1, -1, z),
                coef(1, 0, x), coef(1, 0, y), coef(1, 0, z);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver{SH1D * SH1D.transpose()};
        const auto val_c = solver.eigenvalues();
        auto vec = solver.eigenvectors();
        std::array<double, 3> val;
        for (std::size_t i = 0; i < 3; ++i)
        {
            val[i] = val_c(i) * (vec.col(i).norm());
            vec.col(i).normalize();
        }

        std::array<int, 3> order{0, 1, 2};
        std::sort(order.begin(), order.end(),
                [&val](int i, int j){ return std::abs(val[i]) > std::abs(val[j]); });
        // const Eigen::Vector3d ex = vec.col(order[0]) * (val[order[0]] > 0? 1: -1);
        // const Eigen::Vector3d ez = vec.col(order[2]) * (val[order[2]] > 0? 1: -1);
        // const Eigen::Vector3d ey = ez.cross(ex);
        // Rot.col(0) = ex;
        // Rot.col(1) = ey;
        // Rot.col(2) = ez;
        //
        // sphere.rotate(Rot.transpose());
        Rot(0, 0) = std::abs(val[order[0]]);
        Rot(1, 1) = std::abs(val[order[1]]);
        Rot(2, 2) = std::abs(val[order[2]]);
        Rot = Rot * SH1D.inverse();
std::cout << "SH1D" << SH1D << std::endl;
std::cout << "Rot" << Rot << std::endl;
std::cout << sphere.vertex[0] << std::endl;
        sphere.rotate(Rot);
std::cout << Rot * SH1D << std::endl;
std::cout << sphere.vertex[0] << std::endl;
    }
}; // end of namespace shape


}
