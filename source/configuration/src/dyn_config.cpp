#include "dyn_config.hpp"

using namespace SIMUG;

template <typename RT>
DynConfig<RT>::DynConfig()
{
    this->config_type = conf::type::dyn;
};

template <typename RT>
DynConfig<RT>::DynConfig(const std::string& json_path_):
    Configuration<RT>(json_path_)
{
    this->config_type = conf::type::dyn;
    std::ifstream ifs(this->json_filename);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    nlohmann::json j_dyn = j_input[conf::jsonkey[this->config_type]];

    // Parse all dyn data
    if (!j_dyn["scheme"].empty())
        scheme = dyn::momschnotation[j_dyn["scheme"]];

    if (!j_dyn["pressure"].empty())
        pressure = dyn::pressnotation[j_dyn["pressure"]];

    if (!j_dyn["boundary conditions"].empty())
        boundary_conditions = dyn::bcnotation[j_dyn["boundary conditions"]];
        
    if (!j_dyn["real parameters"].empty())
    {
        std::map<std::string, RT> data = j_dyn["real parameters"].get<std::map<std::string, RT>>();
        for (const auto& [key, val]: data)
            realparameters[key] = val;
    }  

    if (!j_dyn["int parameters"].empty())
    {
        std::map<std::string, int> data = j_dyn["int parameters"].get<std::map<std::string, int>>();
        for (const auto& [key, val]: data)
            intparameters[key] = val;
    }  
};

template <typename RT>
void DynConfig<RT>::SetScheme(dyn::scheme scheme_notation_)
{
    scheme = scheme_notation_;
};

template <typename RT>
void DynConfig<RT>::SetPress(dyn::press press_notation_)
{
    pressure = press_notation_;
};

template <typename RT>
void DynConfig<RT>::SetBC(dyn::bc bc_notation_)
{
    boundary_conditions = bc_notation_;
};

template <typename RT>
void DynConfig<RT>::SetRealParam(const std::string& param_name_,
                                                 RT param_value_)
{
    realparameters[param_name_] = param_value_;
};

template <typename RT>
void DynConfig<RT>::SetIntParam(const std::string& param_name_,
                                               int param_value_)
{
    intparameters[param_name_] = param_value_;
};

template <typename RT>
dyn::scheme DynConfig<RT>::GetScheme() const
{
    return scheme;
};

template <typename RT>
dyn::press DynConfig<RT>::GetPress() const
{
    return pressure;
};

template <typename RT>
inline dyn::bc DynConfig<RT>::GetBC() const
{
    return boundary_conditions;
};

template <typename RT>
RT DynConfig<RT>::GetRealParam(const std::string& param_name_)
{
    return realparameters[param_name_];
};

template <typename RT>
int DynConfig<RT>::GetIntParam(const std::string& param_name_)
{
    return intparameters[param_name_];
};

template <typename RT>
conf::type DynConfig<RT>::GetConfigType() const
{
    return this->config_type;
};

template <typename RT>
void DynConfig<RT>::Log(std::ostream& os) const
{
    os << "===================================================\n";
    os << "MOMENTUM VARIABLES:\n";
    os << "1. scheme: " << dyn::momschname[scheme] << ";\n";
    os << "2. pressure parametrization: " << dyn::pressname[pressure] << ";\n";
    os << "3. boundary conditions: " << dyn::bcname[boundary_conditions] << ";\n";
    int i = 4;
    for (const auto& [key, val]: realparameters)
    {
        os << i << ". " << key << ": " << val << ";\n";
        ++i;
    }
    for (const auto& [key, val]: intparameters)
    {
        os << i << ". " << key << ": " << val << ";\n";
        ++i;
    }
    os << "===================================================\n";
};

// template class instantations
template class DynConfig<float>;
template class DynConfig<double>;