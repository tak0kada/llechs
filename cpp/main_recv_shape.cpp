#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include "cnthd/mesh.hpp"
#include "spherical_harmonics/sph.hpp"
#include "unitsp.hpp"
#include "util.hpp"
#include "read_coef.hpp"

int main(const int argc, const char* const argv[])
{
    //--------------------------------------------------------------------------
    // parse options
    //--------------------------------------------------------------------------
    namespace po = boost::program_options;

    po::options_description opt("");
    opt.add_options()
        ("coef", po::value<std::string>(), ".coef file")
        ("sphere", po::value<std::string>(), ".obj file")
        ("output", po::value<std::string>(), ".obj file");

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

    if (!vm.count("coef") || !vm.count("sphere") || !vm.count("output"))
    {
        std::cerr << opt << std::endl;
        return 1;
    }

    std::string coef_path{};
    std::string sphere_path{};
    std::string output_path{};
    try
    {
        coef_path = vm["coef"].as<std::string>();
        sphere_path  = vm["sphere"].as<std::string>();
        output_path = vm["output"].as<std::string>();
    }
    catch (const po::error_with_option_name& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    if (coef_path.substr(coef_path.find_last_of(".") + 1) != "coef" ||
        sphere_path.substr(sphere_path.find_last_of(".") + 1) != "obj" ||
        output_path.substr(output_path.find_last_of(".") + 1) != "obj")
    {
        std::cerr << opt << std::endl;
        return 1;
    }

    //--------------------------------------------------------------------------
    // recover coordinate
    //--------------------------------------------------------------------------
    using namespace cnthd;
    std::vector<std::vector<double>> coef = read_coef(coef_path);
    Mesh sphere = read_obj(sphere_path);
    to_unit_sphere(sphere);
    Mesh output = sphere;

    for (std::size_t i = 0; i < output.vertex.size(); ++i)
    {
        for (std::size_t ci = 0; ci < 3; ++ci)
        {
            double w{};
            int l{};
            int m{};
            for (const auto c: coef[ci])
            {
                const Vertex v{sphere.vertex[i]};
                w += c * sph::sph_harm(l, m, theta(v), phi(v));
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
            output.vertex[i].p[ci] = w;
        }
    }

    write_obj(output, output_path);

    return 0;
}
