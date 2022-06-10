#include "numerical.hpp"

using namespace SIMUG;

std::vector<double> SIMUG::FromBariToOrdinary(const std::vector<std::vector<double>>& trnodes,
                                       const std::vector<double>& bcoords)
{
    double x0 = trnodes[0][0];
    double y0 = trnodes[0][1];

    double x1 = trnodes[1][0];
    double y1 = trnodes[1][1];

    double x2 = trnodes[2][0];
    double y2 = trnodes[2][1];

    double b0 = bcoords[0];
    double b1 = bcoords[1];
    double b2 = bcoords[2];

    return {x0*b0 + x1*b1 + x2*b2, y0*b0 + y1*b1 + y2*b2};
}

double SIMUG::integral_over_triangle(const std::vector<std::vector<double>>& trnodes,
                              ScalarFunction f)
{
    double w0 = 0.205950504760887;
    double w1 = 0.063691414286223;

    double e00 = 0.124949503233232;
    double e01 = 0.437525248383384;

    double e10 = 0.797112651860071;
    double e11 = 0.165409927389841;
    double e12 = 0.037477420750088;

    double x0 = trnodes[0][0];
    double y0 = trnodes[0][1];

    double x1 = trnodes[1][0];
    double y1 = trnodes[1][1];

    double x2 = trnodes[2][0];
    double y2 = trnodes[2][1];

    double num_int = 0.0;

    // multiplicity = 3
    std::vector<double> p30 = FromBariToOrdinary(trnodes, {e01, e01, e00});
    std::vector<double> p31 = FromBariToOrdinary(trnodes, {e01, e00, e01});
    std::vector<double> p32 = FromBariToOrdinary(trnodes, {e00, e01, e01});

    num_int =(f.Evaluate(p30[0], p30[1]) +
              f.Evaluate(p31[0], p31[1]) + 
              f.Evaluate(p32[0], p32[1]))*w0;

    // multiplicity = 6
    std::vector<double> p60 = FromBariToOrdinary(trnodes, {e10, e11, e12});
    std::vector<double> p61 = FromBariToOrdinary(trnodes, {e10, e12, e11});
    std::vector<double> p62 = FromBariToOrdinary(trnodes, {e11, e10, e12});
    std::vector<double> p63 = FromBariToOrdinary(trnodes, {e11, e12, e10});
    std::vector<double> p64 = FromBariToOrdinary(trnodes, {e12, e10, e11});
    std::vector<double> p65 = FromBariToOrdinary(trnodes, {e12, e11, e10});

    num_int +=(f.Evaluate(p60[0], p60[1]) +
               f.Evaluate(p61[0], p61[1]) + 
               f.Evaluate(p62[0], p62[1]) +
               f.Evaluate(p63[0], p63[1]) +
               f.Evaluate(p64[0], p64[1]) +
               f.Evaluate(p65[0], p65[1]))*w1;
    
    num_int *= SIMUG::trian_square<double>({x0, y0}, {x1, y1}, {x2, y2});
    return num_int;
}