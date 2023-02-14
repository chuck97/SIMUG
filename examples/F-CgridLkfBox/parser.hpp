#pragma once

#include "simug.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

struct Parser
{
    Parser(const std::string& config_file_path);

    double time_step_seconds;
    double total_time_seconds;
    std::string grid_file;
    int output_frequency;
    double alpha_mEVP;
    double beta_mEVP;
    double alpha_stab;
    int Nits_mEVP;
    std::string output_prefix;
    std::string output_dir;
};

Parser::Parser(const std::string& config_file_path)
{
    std::ifstream input_file(config_file_path); 
	
	std::string line;
    std::string str;
	
	// parse  time step
	std::getline(input_file, line);
	std::stringstream ss(line);
	ss >> str;
    if (str != "time_step_seconds")
    {
        SIMUG_ERR("1 line should be: time_step_seconds!");
    }
	ss >> time_step_seconds;
    ss = std::stringstream();

    // parse  total time
	std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "total_time_seconds")
    {
        SIMUG_ERR("2 line should be: total_time_seconds!");
    }
	ss >> total_time_seconds;
    ss = std::stringstream();

    // parse grid path
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "grid_file")
    {
        SIMUG_ERR("3 line should be: grid_file!");
    }
	ss >> grid_file;
    ss = std::stringstream();

    // parse output frequency
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "output_frequency")
    {
        SIMUG_ERR("4 line should be: output_frequency!");
    }
	ss >> output_frequency;
    ss = std::stringstream();

    // parse alpha mEVP
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "alpha_mEVP")
    {
        SIMUG_ERR("5 line should be: alpha_mEVP!");
    }
	ss >> alpha_mEVP;
    ss = std::stringstream();

    // parse beta mEVP
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "beta_mEVP")
    {
        SIMUG_ERR("6 line should be: beta_mEVP!");
    }
	ss >> beta_mEVP;
    ss = std::stringstream();

    // parse alpha stab
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "alpha_stab")
    {
        SIMUG_ERR("7 line should be: alpha_stab!");
    }
	ss >> alpha_stab;
    ss = std::stringstream();

    // parse Nits mEVP
    std::getline(input_file, line);
    ss = std::stringstream(line);
	ss >> str;
    if (str != "Nits_mEVP")
    {
        SIMUG_ERR("8 line should be: Nits_mEVP!");
    }
	ss >> Nits_mEVP;
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