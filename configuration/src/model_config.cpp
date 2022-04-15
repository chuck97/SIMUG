#include "model_config.hpp"

using namespace SIMUG;

template<typename RT>
ModelConfig<RT>::ModelConfig()
{
    this->config_type = conf::type::model;
};

template<typename RT>
ModelConfig<RT>::ModelConfig(const Var<RT>& time_step_, 
                             const Var<RT>& total_time_,
                             const Var<RT>& mesh_)
{
    this->config_type = conf::type::model;
    ModelVars[model::var::dt] = time_step_;
    ModelVars[model::var::T] = total_time_;
    ModelVars[model::var::mesh] = mesh_;
};

template<typename RT>
ModelConfig<RT>::ModelConfig(const Var<RT>& time_step_, 
                             const Var<RT>& total_time_,
                             const Var<RT>& amin_,
                             const Var<RT>& hmin_,
                             const Var<RT>& mesh_)
{
    this->config_type = conf::type::model;
    ModelVars[model::var::dt] = time_step_;
    ModelVars[model::var::T] = total_time_;
    ModelVars[model::var::amin] = amin_;
    ModelVars[model::var::hmin] = hmin_;
    ModelVars[model::var::mesh] = mesh_;
};

template<typename RT>
ModelConfig<RT>::ModelConfig(const std::string& json_path_):
    Configuration<RT>(json_path_)
{
    this->config_type = conf::type::model;
    std::ifstream ifs(this->json_filename);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    nlohmann::json j_model = j_input[conf::jsonkey[conf::type::model]];

    // Parse all model data
    for (model::var item: model::vars)
    {
        if (!j_model[model::name[item]].empty())
        {
            Var<RT> cur_var = {j_model[model::name[item]][0], 
                               j_model[model::name[item]][1]};
            
            ModelVars[item] = cur_var;
        }
    }  
};

template<typename RT>
void ModelConfig<RT>::SetVar(model::var var_notation_,
                             const Var<RT>& var_)
{
    ModelVars[var_notation_] = var_;
};

template<typename RT>
void ModelConfig<RT>::SetVars(const std::map<model::var, Var<RT>>& vars_)
{
    ModelVars = vars_;
};

template<typename RT>
Var<RT> ModelConfig<RT>::GetVar(model::var var_notation_) const
{
    return ModelVars.at(var_notation_);
};

template<typename RT>
conf::type ModelConfig<RT>::GetConfigType() const
{
    return this->config_type;
};

template<typename RT>
void ModelConfig<RT>::Log(std::ostream& os) const
{
    os << "===================================================\n";
    os << "MODEL VARIABLES:\n";
    int i = 1;
    for (const auto& [key, val]: ModelVars)
    {
        os << i << ". " << model::name[key] << " = " << val.value << " " << val.unit << ";\n";
        ++i;
    }
    os << "===================================================\n";
};

// templates instantiation
template class ModelConfig<float>;
template class ModelConfig<double>;