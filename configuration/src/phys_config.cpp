#include "phys_config.hpp"

using namespace SIMUG;

template<typename RT>
PhysConfig<RT>::PhysConfig()
{
    this->config_type = conf::type::phys;
};

template<typename RT>
PhysConfig<RT>::PhysConfig(const std::string& json_path_):
    Configuration<RT>(json_path_)
{
    this->config_type = conf::type::phys;
    std::ifstream ifs(this->json_filename);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    nlohmann::json j_phys = j_input[conf::jsonkey[conf::type::phys]];

    // Parse all phys data
    for (phys::var item: phys::vars)
    {
        if (!j_phys[phys::name[item]].empty())
        {
            Var<RT> cur_var = {j_phys[phys::name[item]][0], 
                               j_phys[phys::name[item]][1]};
            
            PhysVars[item] = cur_var;
        }
    }    
};

template<typename RT>
void PhysConfig<RT>::SetVar(phys::var var_notation_,
                            const Var<RT>& var_)
{
    PhysVars[var_notation_] = var_;
};

template<typename RT>
void PhysConfig<RT>::SetVars(const std::map<phys::var, Var<RT>>& vars_)
{
    PhysVars = vars_;
};

template<typename RT>
Var<RT> PhysConfig<RT>::GetVar(phys::var var_notation_) const
{
    return PhysVars.at(var_notation_);
};

template<typename RT>
conf::type PhysConfig<RT>::GetConfigType() const
{
    return this->config_type;
};

template<typename RT>
void PhysConfig<RT>::Log(std::ostream& os) const
{
    os << "===================================================\n";
    os << "PHYSICAL VARIABLES:\n";
    int i = 1;
    for (const auto& [key, val]: PhysVars)
    {
        os << i << ". " << phys::name[key] << " = " << val.value << " " << val.unit << ";\n";
        ++i;
    }
    os << "===================================================\n";
};

// templates instantiation
template class PhysConfig<float>;
template class PhysConfig<double>;