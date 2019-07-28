#pragma once


#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>

std::vector<double> split(const std::string& str, char delimiter)
{
    std::istringstream iss{str};
    std::vector<double> result;

    std::string buf;
    while (std::getline(iss, buf, delimiter))
    {
        result.push_back(std::stod(buf));
    }
    return result;
}

std::vector<std::vector<double>> read_coef(std::string coef_path)
{
    std::ifstream ifs(coef_path);
    if (!ifs.good())
    {
        std::cerr << "Error can not open " + coef_path << std::endl;
        std::exit(EXIT_FAILURE);
    }


    std::vector<std::vector<double>> coef(3);
    for (int i = 0; i < 3; ++i)
    {
        std::string buf;
        std::getline(ifs, buf);
        coef[i] = split(buf, ' ');
    }

    return coef;
}
