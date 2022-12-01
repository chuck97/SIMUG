#pragma once

#include "vecmath.hpp"
#include "defines.hpp"

namespace SIMUG
{
    enum class FuncType 
    {
         constant,
         linear,
         quadratic
    };

    class ScalarFunction
    {
    public:
        ScalarFunction(FuncType ft_,
                       std::vector<double> params_,
                       std::vector<std::vector<double>> trian_);
        double Evaluate(double x, double y) const;
        FuncType GetType() const;
        std::vector<std::vector<double>> GetTrian() const;
        std::vector<double> GetParams() const;
    private:
        FuncType ft;
        std::vector<double> params;
        std::vector<std::vector<double>> trian;
    };

    class VectorFunction
    {
    public:
        VectorFunction(FuncType ft_,
                       std::vector<double> params1_,
                       std::vector<double> params2_,
                       std::vector<std::vector<double>> trian_);
        FuncType GetType() const;
        std::vector<std::vector<double>> GetTrian() const;
        std::vector<std::vector<double>>  GetParams() const;
    private:
        FuncType ft;
        std::vector<double> params1;
        std::vector<double> params2;
        std::vector<std::vector<double>> trian;
    };

    ScalarFunction operator*(const VectorFunction& v1, const VectorFunction& v2);

    VectorFunction Gradient(ScalarFunction s);

    ScalarFunction Divirgence(VectorFunction v);

    ScalarFunction d_dx(ScalarFunction s);

    ScalarFunction d_dy(ScalarFunction s);

    ScalarFunction operator*(const ScalarFunction& s1, const ScalarFunction& s2);

    VectorFunction operator*(const ScalarFunction& s, const VectorFunction& v);

    ScalarFunction operator+(const ScalarFunction& s1, const ScalarFunction& s2);
}