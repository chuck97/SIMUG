#pragma once

#include "simug.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

enum initScalar
{
    GH, // gaussian hills
    SC  // slotted cylinders
};

enum velField
{
    ND1, // non-div velocity field
    ND2, // non-div velocity field with replacement
    D    // div velocty field
};

struct Parser
{
    Parser(const std::string& config_file_path);
    double courant_number;
    std::string grid_file;
    int output_frequency;
    initScalar init_scalar;
    velField  vel_field;
    SIMUG::adv::timeScheme time_scheme;
    SIMUG::adv::spaceScheme space_scheme;
    SIMUG::adv::advFilter adv_filter;
    std::string output_prefix;
    std::string output_dir;
};

initScalar GetInitScal(const std::string& name)
{
    if (name == (std::string)"gaussian_hills")
    {
        return initScalar::GH;
    }
    else if (name == (std::string)"slotted_cylinders")
    {
        return initScalar::SC;
    }
    else
    {
        SIMUG_ERR("possible variants for initial scalar distribution: gaussian_hills, slotted_cylinders!");
    }
}

velField GetVelField(const std::string& name)
{
    if (name == (std::string)"non_div_1")
    {
        return velField::ND1;
    }
    else if (name == (std::string)"non_div_2")
    {
        return velField::ND2;
    }
    else if (name == (std::string)"div")
    {
        return velField::D;
    }
    else
    {
        SIMUG_ERR("possible variants for velocity field: non_div_1, non_div_2, div!");
    }
}

SIMUG::adv::timeScheme GetTimeScheme(const std::string& name)
{
    if (name == (std::string)"Euler")
    {
        return SIMUG::adv::timeScheme::Euler;
    }
    else if (name == (std::string)"TRK2")
    {
        return SIMUG::adv::timeScheme::TRK2;
    }
    else
    {
        SIMUG_ERR("possible variants for advection time scheme on C-grid: Euler, TRK2!");
    }
}

SIMUG::adv::spaceScheme GetSpaceScheme(const std::string& name)
{
    if (name == (std::string)"upwind")
    {
        return SIMUG::adv::spaceScheme::FVupwind;
    }
    else if (name == (std::string)"MUST")
    {
        return SIMUG::adv::spaceScheme::MUST;
    }
    else if (name == (std::string)"MUSCL")
    {
        return SIMUG::adv::spaceScheme::MUSCL;
    }
    else
    {
        SIMUG_ERR("possible variants for advection space scheme on C-grid: upwind, MUST, MUSCL!");
    }
}

SIMUG::adv::advFilter GetFilter(const std::string& name)
{
    if (name == (std::string)"none")
    {
        return SIMUG::adv::advFilter::none;
    }
    else if (name == (std::string)"minmod")
    {
        return SIMUG::adv::advFilter::Minmod;
    }
    else if (name == (std::string)"VanLeer")
    {
        return SIMUG::adv::advFilter::VanLeer;
    }
    else if (name == (std::string)"Superbee")
    {
        return SIMUG::adv::advFilter::Superbee;
    }
    else if (name == (std::string)"BarthJesperson")
    {
        return SIMUG::adv::advFilter::BarthJesperson;
    }
    else
    {
        SIMUG_ERR("possible variants for advection filter on C-grid: none, minmod, VanLeer, Superbee, BarthJesperson!");
    }
}

Parser::Parser(const std::string& config_file_path)
{
    std::ifstream input_file(config_file_path); 
	
	std::string line;
    std::string str;
	
	// parse Courant number
	std::getline(input_file, line);
	std::stringstream ss(line);
	ss >> str;
    if (str != "courant_number")
    {
        SIMUG_ERR("1 line should be: courant_number!");
    }
	ss >> courant_number;
    ss = std::stringstream();

    // parse grid path
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "grid_file")
    {
        SIMUG_ERR("2 line should be: grid_file!");
    }
	ss >> grid_file;
    ss = std::stringstream();

    // parse output frequency
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "output_frequency")
    {
        SIMUG_ERR("3 line should be: output_frequency!");
    }
	ss >> output_frequency;
    ss = std::stringstream();

    // parse init scalar
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "initial_scalar")
    {
        SIMUG_ERR("4 line should be: initial_scalar!");
    }
	ss >> str;
    init_scalar = GetInitScal(str);
    ss = std::stringstream();

    // parse vel field
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "velocity_field")
    {
        SIMUG_ERR("5 line should be: velocity_field!");
    }
	ss >> str;
    vel_field = GetVelField(str);
    ss = std::stringstream();

    // parse advecton time scheme
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "time_scheme")
    {
        SIMUG_ERR("6 line should be: time_scheme!");
    }
	ss >> str;
    time_scheme = GetTimeScheme(str);
    ss = std::stringstream();

    // parse advecton space scheme
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "space_scheme")
    {
        SIMUG_ERR("7 line should be: space_scheme!");
    }
	ss >> str;
    space_scheme = GetSpaceScheme(str);
    ss = std::stringstream();

    // parse advecton filter
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "adv_filter")
    {
        SIMUG_ERR("8 line should be: adv_filter!");
    }
	ss >> str;
    adv_filter = GetFilter(str);
    ss = std::stringstream();

    // parse output prefix
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "output_prefix")
    {
        SIMUG_ERR("9 line should be: output_prefix!");
    }
	ss >> output_prefix;
    ss = std::stringstream();

    // parse output dir
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "output_dir")
    {
        SIMUG_ERR("10 line should be: output_dir!");
    }
	ss >> output_dir;
};