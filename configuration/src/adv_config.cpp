#include "adv_config.hpp"

using namespace SIMUG;

template <typename RT>
AdvConfig<RT>::AdvConfig()
{
    this->config_type = conf::type::adv;
};

template <typename RT>
AdvConfig<RT>::AdvConfig(const std::string& json_path_):
    Configuration<RT>(json_path_)
{
    this->config_type = conf::type::adv;
    std::ifstream ifs(this->json_filename);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    nlohmann::json j_adv = j_input[conf::jsonkey[conf::type::adv]];

    // Parse all advection data
    if (!j_adv["scheme"].empty())
        scheme = adv::advschnotation[j_adv["scheme"]];

    if (!j_adv["filter"].empty())
        filter = adv::advfiltnotation[j_adv["filter"]];

    if (!j_adv["real parameters"].empty())
    {
        std::map<std::string, RT> data = j_adv["real parameters"].get<std::map<std::string, RT>>();
        for (const auto& [key, val]: data)
            realparameters[key] = val;
    }

    if (!j_adv["int parameters"].empty())
    {
        std::map<std::string, int> data = j_adv["int parameters"].get<std::map<std::string, int>>();
        for (const auto& [key, val]: data)
            intparameters[key] = val;
    }  
};

template <typename RT>
void AdvConfig<RT>::SetScheme(adv::scheme scheme_notation_)
{
    scheme = scheme_notation_;
};

template <typename RT>
void AdvConfig<RT>::SetFilter(adv::filter filter_notation_)
{
    filter = filter_notation_;
};

template <typename RT>
void AdvConfig<RT>::SetRealParam(const std::string& param_name_,
                                                RT  param_value_)
{
    realparameters[param_name_] = param_value_;
};

template <typename RT>
void AdvConfig<RT>::SetIntParam(const std::string& param_name_,
                                              int  param_value_)
{
    intparameters[param_name_] = param_value_;
};

template <typename RT>
adv::scheme AdvConfig<RT>::GetScheme() const
{
    return scheme;
};

template <typename RT>
adv::filter AdvConfig<RT>::GetFilter() const
{
    return filter;
};

template <typename RT>
RT AdvConfig<RT>::GetRealParam(const std::string param_name_)
{
    return realparameters[param_name_];
};

template <typename RT>
int AdvConfig<RT>::GetIntParam(const std::string param_name_)
{
    return intparameters[param_name_];
};

template <typename RT>
conf::type AdvConfig<RT>::GetConfigType() const
{
    return this->config_type;
};

template <typename RT>
void AdvConfig<RT>::Log(std::ostream& os) const
{
    os << "===================================================\n";
    os << "ADVECTION VARIABLES:\n";
    os << "1. scheme: " << adv::advschname[scheme] << ";\n";
    os << "2. filter: " << adv::advfiltname[filter] << ";\n";
    int i = 3;
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
template class AdvConfig<float>;
template class AdvConfig<double>;