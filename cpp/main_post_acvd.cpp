#include <utility>
#include <boost/program_options.hpp>
#include "cnthd/mesh.hpp"
#include "post_acvd.hpp"

int main(int argc, char *argv[])
{
    //--------------------------------------------------------------------------
    // check command line options
    //--------------------------------------------------------------------------
    std::string in_n{};
    std::string out_n{};

    {
        namespace po = boost::program_options;

        po::options_description opt("");
        opt.add_options()
            ("input", po::value<std::string>(), ".ply: input file")
            ("output", po::value<std::string>(), ".obj: output file");

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, opt), vm);
        }
        catch (const po::error_with_option_name& e)
        {
            std::cout << e.what() << std::endl;
            return 1;
        }
        vm.notify();

        if (!vm.count("input") || !vm.count("output"))
        {
            std::cerr << opt << std::endl;
            return 1;
        }

        try
        {
            in_n = vm["input"].as<std::string>();
            out_n = vm["output"].as<std::string>();
        }
        catch (const po::error_with_option_name& e)
        {
            std::cout << e.what() << std::endl;
            return 1;
        }

        // file extensions: they must be ".obj" files
        std::string in_ext{
            in_n.substr(in_n.find_last_of(".") + 1)};
        std::string output_ext{
            out_n.substr(out_n.find_last_of(".") + 1)};

        bool fail{false};
        if (in_ext != "ply")
        {
            std::cerr << "ERROR: Input file must be an ply file." << std::endl;
            fail = true;
        }
        if (output_ext != "obj") {
            std::cerr << "ERROR: Output file must be an obj file." << std::endl;
            fail = true;
        }
        if (fail)
        {
            return 1;
        }

    } // end parse

    //--------------------------------------------------------------------------
    // check command line options
    //--------------------------------------------------------------------------
    const auto [array, face] = post_acvd::read_ply(in_n);
    cnthd::Mesh mesh{array, face};
    mesh.fix_orientation();
    cnthd::write_obj(mesh, out_n);

    return 0;
}
