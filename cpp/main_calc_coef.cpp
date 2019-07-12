#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <boost/program_options.hpp>
#include "cnthd/mesh.hpp"
#include "spherical_harmonics/sph.hpp"
#include "calc_coef.hpp"
#include "unitsp.hpp"

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    if (vec.empty())
    {
        os << "";
    }
    else{
        for (std::size_t i = 0; i < vec.size() - 1; ++i)
        {
            os << vec[i] << " ";
        }
        os << vec[vec.size() - 1];
    }
    return os;
}

int main(int argc, char *argv[])
{
    //--------------------------------------------------------------------------
    // check command line options
    //--------------------------------------------------------------------------
    std::string xyz_n{};
    // NOTE: input sphere must have the same ordering of vertices as in xyz
    std::string sphere_n{};
    std::string output_n{};
    unsigned int L{};

    {
        namespace po = boost::program_options;

        po::options_description opt("");
        opt.add_options()
            ("original", po::value<std::string>(), ".obj: original file")
            ("sphere", po::value<std::string>(), ".obj: file mapped from the original to sphere")
            ("output", po::value<std::string>(), ".coef: output coef file")
            ("sph_degree", po::value<unsigned int>(), "int: maximum L of sphere harmonics");

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, opt), vm);
        }
        catch (const po::error_with_option_name& e)
        {
            std::cout << e.what() << std::endl;
        }
        vm.notify();

        if (!vm.count("original") || !vm.count("sphere") || !vm.count("output") || !vm.count("sph_degree"))
        {
            std::cerr << opt << std::endl;
            return 1;
        }

        try
        {
            xyz_n = vm["original"].as<std::string>();
            sphere_n = vm["sphere"].as<std::string>();
            output_n = vm["output"].as<std::string>();
            L = vm["sph_degree"].as<unsigned int>();
        }
        catch (const po::error_with_option_name& e)
        {
            std::cout << e.what() << std::endl;
            return 1;
        }

        // file extensions: they must be ".obj" files
        std::string xyz_ext{
            xyz_n.substr(xyz_n.find_last_of(".") + 1)};
        std::string sphere_ext{
            sphere_n.substr(sphere_n.find_last_of(".") + 1)};
        std::string output_ext{
            output_n.substr(output_n.find_last_of(".") + 1)};

        bool fail{false};
        if (xyz_ext != "obj" || sphere_ext != "obj")
        {
            std::cerr << "ERROR: Input file must be an obj file." << std::endl;
            fail = true;
        }
        if (output_ext != "coef") {
            std::cerr << "ERROR: Output file must be a coef file." << std::endl;
            fail = true;
        }
        if (fail)
        {
            return 1;
        }

    } // end parse

    //--------------------------------------------------------------------------
    // main program below
    //--------------------------------------------------------------------------
    cnthd::Mesh xyz = cnthd::read_obj(xyz_n);
    cnthd::Mesh sphere = cnthd::read_obj(sphere_n);
    assert(xyz.nE == sphere.nE);
    assert(xyz.nF == sphere.nF);
    assert(xyz.nV == sphere.nV);
    to_unit_sphere(sphere);

    if (xyz.num_genus() != 0)
    {
        std::cerr << "ERROR: Genus of this mesh is not zero." << std::endl;
        return 1;
    }

    std::ofstream ofs(output_n);
    if (!ofs.good())
    {
        std::cerr << "Error: Cannot open coef file!" << std::endl;
        return 1;
    }

    cnthd::shape::shmap M{xyz, sphere};
    // M.align_domain();
    const std::array<std::vector<double>, 3> coef{M.calc_coef((L + 1)*(L + 1))};
    for (std::size_t i = 0; i < 3; ++i)
    {
        ofs << coef[i] << std::endl;
    }
    ofs.close();

    return 0;
}
