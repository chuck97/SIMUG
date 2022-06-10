#include "local_assembling.hpp"

using namespace SIMUG;

std::vector<std::vector<double>> SIMUG::LocalMassMatrixAssembling(const std::vector<std::vector<double>>& node_coords)
{
    double x0 = node_coords[0][0];
    double y0 = node_coords[0][1];

    double x1 = node_coords[1][0];
    double y1 = node_coords[1][1];

    double x2 = node_coords[2][0];
    double y2 = node_coords[2][1];

    std::vector<double> coeffs0 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{1.0, 0.0, 0.0});
    
    std::vector<double> coeffs1 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 1.0, 0.0});

    std::vector<double> coeffs2 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 0.0, 1.0});

    ScalarFunction phi0(FuncType::linear, coeffs0, node_coords);
    ScalarFunction phi1(FuncType::linear, coeffs1, node_coords);
    ScalarFunction phi2(FuncType::linear, coeffs2, node_coords);

    double phi0_phi0 = integral_over_triangle(node_coords, phi0*phi0);
    double phi0_phi1 = integral_over_triangle(node_coords, phi0*phi1);
    double phi0_phi2 = integral_over_triangle(node_coords, phi0*phi2);

    double phi1_phi0 = integral_over_triangle(node_coords, phi1*phi0);
    double phi1_phi1 = integral_over_triangle(node_coords, phi1*phi1);
    double phi1_phi2 = integral_over_triangle(node_coords, phi1*phi2);

    double phi2_phi0 = integral_over_triangle(node_coords, phi2*phi0);
    double phi2_phi1 = integral_over_triangle(node_coords, phi2*phi1);
    double phi2_phi2 = integral_over_triangle(node_coords, phi2*phi2);

    return {{phi0_phi0, phi0_phi1, phi0_phi2},
            {phi1_phi0, phi1_phi1, phi1_phi2},
            {phi2_phi0, phi2_phi1, phi2_phi2}};
}

double SIMUG::L_ij_entry(const std::vector<std::vector<double>>& node_coords,
                      const ScalarFunction& phi_i,
                      const ScalarFunction& phi_j,
                      const ScalarFunction& u0,
                      const ScalarFunction& u1)
{
    double term1, term2, term3;

    // first set
    term1 = integral_over_triangle(node_coords, d_dx(phi_i)*d_dx(phi_j)*u0*u0) +
        2.0*integral_over_triangle(node_coords, d_dx(phi_i)*d_dx(u0)*u0*phi_j);
    
    // second set
    term2 = integral_over_triangle(node_coords, d_dy(phi_i)*d_dy(phi_j)*u1*u1) +
        2.0*integral_over_triangle(node_coords, d_dy(phi_i)*d_dy(u1)*u1*phi_j);

    // third set
    term3 = 2.0*integral_over_triangle(node_coords, d_dy(phi_i)*d_dx(phi_j)*u0*u1) +
            2.0*integral_over_triangle(node_coords, d_dy(phi_i)*d_dx(u1)*u0*phi_j) +
            2.0*integral_over_triangle(node_coords, d_dy(phi_i)*d_dx(u0)*u1*phi_j);
    
    return (term1 + term2 + term3);
}

std::vector<double> SIMUG::LocalTG2RhsAssembling(const std::vector<std::vector<double>>& node_coords,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    double time_step)
{
    double u00 =  uvalues[0][0];
    double u01 =  uvalues[0][1];

    double u10 =  uvalues[1][0];
    double u11 =  uvalues[1][1];

    double u20 =  uvalues[2][0];
    double u21 =  uvalues[2][1];

    double x0 = node_coords[0][0];
    double y0 = node_coords[0][1];

    double x1 = node_coords[1][0];
    double y1 = node_coords[1][1];

    double x2 = node_coords[2][0];
    double y2 = node_coords[2][1];

    std::vector<double> coeffs0 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{1.0, 0.0, 0.0});
    
    std::vector<double> coeffs1 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 1.0, 0.0});

    std::vector<double> coeffs2 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 0.0, 1.0});

    ScalarFunction phi0(FuncType::linear, coeffs0, node_coords);
    ScalarFunction phi1(FuncType::linear, coeffs1, node_coords);
    ScalarFunction phi2(FuncType::linear, coeffs2, node_coords);

    double phi0_phi0 = integral_over_triangle(node_coords, phi0*phi0);
    double phi0_phi1 = integral_over_triangle(node_coords, phi0*phi1);
    double phi0_phi2 = integral_over_triangle(node_coords, phi0*phi2);

    double phi1_phi0 = integral_over_triangle(node_coords, phi1*phi0);
    double phi1_phi1 = integral_over_triangle(node_coords, phi1*phi1);
    double phi1_phi2 = integral_over_triangle(node_coords, phi1*phi2);

    double phi2_phi0 = integral_over_triangle(node_coords, phi2*phi0);
    double phi2_phi1 = integral_over_triangle(node_coords, phi2*phi1);
    double phi2_phi2 = integral_over_triangle(node_coords, phi2*phi2);

    std::vector<std::vector<double>> M = {{phi0_phi0, phi0_phi1, phi0_phi2},
                                          {phi1_phi0, phi1_phi1, phi1_phi2},
                                          {phi2_phi0, phi2_phi1, phi2_phi2}};
    // make u vector function
    ScalarFunction u00_comp(FuncType::constant, {u00}, node_coords);
    ScalarFunction u01_comp(FuncType::constant, {u01}, node_coords);
    ScalarFunction u10_comp(FuncType::constant, {u10}, node_coords);
    ScalarFunction u11_comp(FuncType::constant, {u11}, node_coords);
    ScalarFunction u20_comp(FuncType::constant, {u20}, node_coords);
    ScalarFunction u21_comp(FuncType::constant, {u21}, node_coords);

    ScalarFunction u0 = u00_comp*phi0 + u10_comp*phi1 + u20_comp*phi2;
    ScalarFunction u1 = u01_comp*phi0 + u11_comp*phi1 + u21_comp*phi2;

    VectorFunction u(FuncType::linear, u0.GetParams(), u1.GetParams(), node_coords);

    // K matrix assembling

    double K00, K01, K02, K10, K11, K12, K20, K21, K22; 
    
    K00 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi0));
    K01 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi0));
    K02 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi0));
    K10 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi1));
    K11 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi1));
    K12 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi1));
    K20 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi2));
    K21 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi2));
    K22 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi2));

    std::vector<std::vector<double>> K = {{K00, K01, K02},
                                          {K10, K11, K12},
                                          {K20, K21, K22}};

    // N matrix assembling

    double N00, N01, N02, N10, N11, N12, N20, N21, N22; 

    N00 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi0));
    N01 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi0));
    N02 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi0));
    N10 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi1));
    N11 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi1));
    N12 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi1));
    N20 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi2));
    N21 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi2));
    N22 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi2));

    std::vector<std::vector<double>> N = {{N00, N01, N02},
                                          {N10, N11, N12},
                                          {N20, N21, N22}}; 

    std::vector<std::vector<double>> RhsMatrix = M + time_step*K - ((time_step*time_step)/2.0)*N;

    return {RhsMatrix*localmass};
}

std::vector<double> SIMUG::LocalCG2RhsAssembling(const std::vector<std::vector<double>>& node_coords,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    double time_step)
{
    double u00 =  uvalues[0][0];
    double u01 =  uvalues[0][1];

    double u10 =  uvalues[1][0];
    double u11 =  uvalues[1][1];

    double u20 =  uvalues[2][0];
    double u21 =  uvalues[2][1];

    double x0 = node_coords[0][0];
    double y0 = node_coords[0][1];

    double x1 = node_coords[1][0];
    double y1 = node_coords[1][1];

    double x2 = node_coords[2][0];
    double y2 = node_coords[2][1];

    std::vector<double> coeffs0 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{1.0, 0.0, 0.0});
    
    std::vector<double> coeffs1 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 1.0, 0.0});

    std::vector<double> coeffs2 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 0.0, 1.0});

    ScalarFunction phi0(FuncType::linear, coeffs0, node_coords);
    ScalarFunction phi1(FuncType::linear, coeffs1, node_coords);
    ScalarFunction phi2(FuncType::linear, coeffs2, node_coords);

    double phi0_phi0 = integral_over_triangle(node_coords, phi0*phi0);
    double phi0_phi1 = integral_over_triangle(node_coords, phi0*phi1);
    double phi0_phi2 = integral_over_triangle(node_coords, phi0*phi2);

    double phi1_phi0 = integral_over_triangle(node_coords, phi1*phi0);
    double phi1_phi1 = integral_over_triangle(node_coords, phi1*phi1);
    double phi1_phi2 = integral_over_triangle(node_coords, phi1*phi2);

    double phi2_phi0 = integral_over_triangle(node_coords, phi2*phi0);
    double phi2_phi1 = integral_over_triangle(node_coords, phi2*phi1);
    double phi2_phi2 = integral_over_triangle(node_coords, phi2*phi2);

    std::vector<std::vector<double>> M = {{phi0_phi0, phi0_phi1, phi0_phi2},
                                          {phi1_phi0, phi1_phi1, phi1_phi2},
                                          {phi2_phi0, phi2_phi1, phi2_phi2}};
    // make u vector function
    ScalarFunction u00_comp(FuncType::constant, {u00}, node_coords);
    ScalarFunction u01_comp(FuncType::constant, {u01}, node_coords);
    ScalarFunction u10_comp(FuncType::constant, {u10}, node_coords);
    ScalarFunction u11_comp(FuncType::constant, {u11}, node_coords);
    ScalarFunction u20_comp(FuncType::constant, {u20}, node_coords);
    ScalarFunction u21_comp(FuncType::constant, {u21}, node_coords);

    ScalarFunction u0 = u00_comp*phi0 + u10_comp*phi1 + u20_comp*phi2;
    ScalarFunction u1 = u01_comp*phi0 + u11_comp*phi1 + u21_comp*phi2;

    VectorFunction u(FuncType::linear, u0.GetParams(), u1.GetParams(), node_coords);

    // K matrix assembling
    double K00, K01, K02, K10, K11, K12, K20, K21, K22;
     
    K00 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi0));
    K01 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi0));
    K02 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi0));
    K10 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi1));
    K11 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi1));
    K12 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi1));
    K20 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi2));
    K21 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi2));
    K22 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi2));

    std::vector<std::vector<double>> K = {{K00, K01, K02},
                                          {K10, K11, K12},
                                          {K20, K21, K22}};

    // oN matrix assembling
    double N00, N01, N02, N10, N11, N12, N20, N21, N22;  

    N00 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi0));
    N01 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi0));
    N02 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi0));
    N10 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi1));
    N11 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi1));
    N12 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi1));
    N20 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi2));
    N21 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi2));
    N22 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi2));

    std::vector<std::vector<double>> N = {{N00, N01, N02},
                                           {N10, N11, N12},
                                           {N20, N21, N22}}; 
    
    
    std::vector<std::vector<double>> L;

    // L matrix assembling
    double L00 = L_ij_entry(node_coords, phi0, phi0, u0, u1);
    double L01 = L_ij_entry(node_coords, phi1, phi0, u0, u1);
    double L02 = L_ij_entry(node_coords, phi2, phi0, u0, u1);
    double L10 = L_ij_entry(node_coords, phi0, phi1, u0, u1);
    double L11 = L_ij_entry(node_coords, phi1, phi1, u0, u1);
    double L12 = L_ij_entry(node_coords, phi2, phi1, u0, u1);
    double L20 = L_ij_entry(node_coords, phi0, phi2, u0, u1);
    double L21 = L_ij_entry(node_coords, phi1, phi2, u0, u1);
    double L22 = L_ij_entry(node_coords, phi2, phi2, u0, u1);

    L = {{L00, L01, L02},
         {L10, L11, L12},
         {L20, L21, L22}}; 

    std::vector<std::vector<double>> RhsMatrix;
    RhsMatrix = M + time_step*K - ((time_step*time_step)/2.0)*N - ((time_step*time_step)/2.0)*L;
    return {RhsMatrix*localmass};
}

std::vector<double> SIMUG::LocalTTG2RhsAssembling(const std::vector<std::vector<double>>& node_coords,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    const std::vector<double>& localmass_half,
                                                    double time_step,
                                                    int step)
{
    double u00 =  uvalues[0][0];
    double u01 =  uvalues[0][1];

    double u10 =  uvalues[1][0];
    double u11 =  uvalues[1][1];

    double u20 =  uvalues[2][0];
    double u21 =  uvalues[2][1];

    double x0 = node_coords[0][0];
    double y0 = node_coords[0][1];

    double x1 = node_coords[1][0];
    double y1 = node_coords[1][1];

    double x2 = node_coords[2][0];
    double y2 = node_coords[2][1];

    std::vector<double> coeffs0 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{1.0, 0.0, 0.0});
    
    std::vector<double> coeffs1 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 1.0, 0.0});

    std::vector<double> coeffs2 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 0.0, 1.0});

    ScalarFunction phi0(FuncType::linear, coeffs0, node_coords);
    ScalarFunction phi1(FuncType::linear, coeffs1, node_coords);
    ScalarFunction phi2(FuncType::linear, coeffs2, node_coords);

    double phi0_phi0 = integral_over_triangle(node_coords, phi0*phi0);
    double phi0_phi1 = integral_over_triangle(node_coords, phi0*phi1);
    double phi0_phi2 = integral_over_triangle(node_coords, phi0*phi2);

    double phi1_phi0 = integral_over_triangle(node_coords, phi1*phi0);
    double phi1_phi1 = integral_over_triangle(node_coords, phi1*phi1);
    double phi1_phi2 = integral_over_triangle(node_coords, phi1*phi2);

    double phi2_phi0 = integral_over_triangle(node_coords, phi2*phi0);
    double phi2_phi1 = integral_over_triangle(node_coords, phi2*phi1);
    double phi2_phi2 = integral_over_triangle(node_coords, phi2*phi2);

    std::vector<std::vector<double>> M = {{phi0_phi0, phi0_phi1, phi0_phi2},
                                          {phi1_phi0, phi1_phi1, phi1_phi2},
                                          {phi2_phi0, phi2_phi1, phi2_phi2}};
    // make u vector function
    ScalarFunction u00_comp(FuncType::constant, {u00}, node_coords);
    ScalarFunction u01_comp(FuncType::constant, {u01}, node_coords);
    ScalarFunction u10_comp(FuncType::constant, {u10}, node_coords);
    ScalarFunction u11_comp(FuncType::constant, {u11}, node_coords);
    ScalarFunction u20_comp(FuncType::constant, {u20}, node_coords);
    ScalarFunction u21_comp(FuncType::constant, {u21}, node_coords);

    ScalarFunction u0 = u00_comp*phi0 + u10_comp*phi1 + u20_comp*phi2;
    ScalarFunction u1 = u01_comp*phi0 + u11_comp*phi1 + u21_comp*phi2;

    VectorFunction u(FuncType::linear, u0.GetParams(), u1.GetParams(), node_coords);

    // K matrix assembling

    double K00, K01, K02, K10, K11, K12, K20, K21, K22; 
    
    K00 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi0));
    K01 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi0));
    K02 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi0));
    K10 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi1));
    K11 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi1));
    K12 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi1));
    K20 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi2));
    K21 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi2));
    K22 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi2));

    std::vector<std::vector<double>> K = {{K00, K01, K02},
                                          {K10, K11, K12},
                                          {K20, K21, K22}};

    if (step == 1)
    {
        return {M*localmass + (time_step/2.0)*K*localmass};
    }
    else if (step == 2)
    {
        return {M*localmass + time_step*K*localmass_half};
    }
    else
    {
        SIMUG_ERR("only 2 step TG method available");
    }
}


std::vector<double> SIMUG::LocalTTG3RhsAssembling(const std::vector<std::vector<double>>& node_coords,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    const std::vector<double>& localmass_half,
                                                    double time_step,
                                                    int step)
{
    double u00 =  uvalues[0][0];
    double u01 =  uvalues[0][1];

    double u10 =  uvalues[1][0];
    double u11 =  uvalues[1][1];

    double u20 =  uvalues[2][0];
    double u21 =  uvalues[2][1];

    double x0 = node_coords[0][0];
    double y0 = node_coords[0][1];

    double x1 = node_coords[1][0];
    double y1 = node_coords[1][1];

    double x2 = node_coords[2][0];
    double y2 = node_coords[2][1];

    std::vector<double> coeffs0 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{1.0, 0.0, 0.0});
    
    std::vector<double> coeffs1 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 1.0, 0.0});

    std::vector<double> coeffs2 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 0.0, 1.0});

    ScalarFunction phi0(FuncType::linear, coeffs0, node_coords);
    ScalarFunction phi1(FuncType::linear, coeffs1, node_coords);
    ScalarFunction phi2(FuncType::linear, coeffs2, node_coords);

    double phi0_phi0 = integral_over_triangle(node_coords, phi0*phi0);
    double phi0_phi1 = integral_over_triangle(node_coords, phi0*phi1);
    double phi0_phi2 = integral_over_triangle(node_coords, phi0*phi2);

    double phi1_phi0 = integral_over_triangle(node_coords, phi1*phi0);
    double phi1_phi1 = integral_over_triangle(node_coords, phi1*phi1);
    double phi1_phi2 = integral_over_triangle(node_coords, phi1*phi2);

    double phi2_phi0 = integral_over_triangle(node_coords, phi2*phi0);
    double phi2_phi1 = integral_over_triangle(node_coords, phi2*phi1);
    double phi2_phi2 = integral_over_triangle(node_coords, phi2*phi2);

    std::vector<std::vector<double>> M = {{phi0_phi0, phi0_phi1, phi0_phi2},
                                          {phi1_phi0, phi1_phi1, phi1_phi2},
                                          {phi2_phi0, phi2_phi1, phi2_phi2}};
    // make u vector function
    ScalarFunction u00_comp(FuncType::constant, {u00}, node_coords);
    ScalarFunction u01_comp(FuncType::constant, {u01}, node_coords);
    ScalarFunction u10_comp(FuncType::constant, {u10}, node_coords);
    ScalarFunction u11_comp(FuncType::constant, {u11}, node_coords);
    ScalarFunction u20_comp(FuncType::constant, {u20}, node_coords);
    ScalarFunction u21_comp(FuncType::constant, {u21}, node_coords);

    ScalarFunction u0 = u00_comp*phi0 + u10_comp*phi1 + u20_comp*phi2;
    ScalarFunction u1 = u01_comp*phi0 + u11_comp*phi1 + u21_comp*phi2;

    VectorFunction u(FuncType::linear, u0.GetParams(), u1.GetParams(), node_coords);

    // K matrix assembling

    double K00, K01, K02, K10, K11, K12, K20, K21, K22; 
    
    K00 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi0));
    K01 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi0));
    K02 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi0));
    K10 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi1));
    K11 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi1));
    K12 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi1));
    K20 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi2));
    K21 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi2));
    K22 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi2));

    std::vector<std::vector<double>> K = {{K00, K01, K02},
                                          {K10, K11, K12},
                                          {K20, K21, K22}};

    // N matrix assembling

    double N00, N01, N02, N10, N11, N12, N20, N21, N22; 

    N00 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi0));
    N01 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi0));
    N02 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi0));
    N10 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi1));
    N11 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi1));
    N12 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi1));
    N20 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi2));
    N21 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi2));
    N22 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi2));

    std::vector<std::vector<double>> N = {{N00, N01, N02},
                                          {N10, N11, N12},
                                          {N20, N21, N22}}; 

    if (step == 1)
    {
        return{M*localmass + (time_step/3.0)*K*localmass - ((time_step*time_step)/9.0)*N*localmass};
    }
    else
    {
        return{M*localmass + time_step*K*localmass - ((time_step*time_step)/2.0)*N*localmass_half};
    }
}

std::vector<double> SIMUG::LocalTTG4RhsAssembling(const std::vector<std::vector<double>>& node_coords,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    const std::vector<double>& localmass_half,
                                                    double time_step,
                                                    int step)
{
    double u00 =  uvalues[0][0];
    double u01 =  uvalues[0][1];

    double u10 =  uvalues[1][0];
    double u11 =  uvalues[1][1];

    double u20 =  uvalues[2][0];
    double u21 =  uvalues[2][1];

    double x0 = node_coords[0][0];
    double y0 = node_coords[0][1];

    double x1 = node_coords[1][0];
    double y1 = node_coords[1][1];

    double x2 = node_coords[2][0];
    double y2 = node_coords[2][1];

    std::vector<double> coeffs0 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{1.0, 0.0, 0.0});
    
    std::vector<double> coeffs1 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 1.0, 0.0});

    std::vector<double> coeffs2 = solve_linear_system({{x0, y0, 1.0},
                                                       {x1, y1, 1.0},
                                                       {x2, y2, 1.0}},
                                                      std::vector<double>{0.0, 0.0, 1.0});

    ScalarFunction phi0(FuncType::linear, coeffs0, node_coords);
    ScalarFunction phi1(FuncType::linear, coeffs1, node_coords);
    ScalarFunction phi2(FuncType::linear, coeffs2, node_coords);

    double phi0_phi0 = integral_over_triangle(node_coords, phi0*phi0);
    double phi0_phi1 = integral_over_triangle(node_coords, phi0*phi1);
    double phi0_phi2 = integral_over_triangle(node_coords, phi0*phi2);

    double phi1_phi0 = integral_over_triangle(node_coords, phi1*phi0);
    double phi1_phi1 = integral_over_triangle(node_coords, phi1*phi1);
    double phi1_phi2 = integral_over_triangle(node_coords, phi1*phi2);

    double phi2_phi0 = integral_over_triangle(node_coords, phi2*phi0);
    double phi2_phi1 = integral_over_triangle(node_coords, phi2*phi1);
    double phi2_phi2 = integral_over_triangle(node_coords, phi2*phi2);

    std::vector<std::vector<double>> M = {{phi0_phi0, phi0_phi1, phi0_phi2},
                                          {phi1_phi0, phi1_phi1, phi1_phi2},
                                          {phi2_phi0, phi2_phi1, phi2_phi2}};
    // make u vector function
    ScalarFunction u00_comp(FuncType::constant, {u00}, node_coords);
    ScalarFunction u01_comp(FuncType::constant, {u01}, node_coords);
    ScalarFunction u10_comp(FuncType::constant, {u10}, node_coords);
    ScalarFunction u11_comp(FuncType::constant, {u11}, node_coords);
    ScalarFunction u20_comp(FuncType::constant, {u20}, node_coords);
    ScalarFunction u21_comp(FuncType::constant, {u21}, node_coords);

    ScalarFunction u0 = u00_comp*phi0 + u10_comp*phi1 + u20_comp*phi2;
    ScalarFunction u1 = u01_comp*phi0 + u11_comp*phi1 + u21_comp*phi2;

    VectorFunction u(FuncType::linear, u0.GetParams(), u1.GetParams(), node_coords);

    // K matrix assembling

    double K00, K01, K02, K10, K11, K12, K20, K21, K22; 
    
    K00 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi0));
    K01 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi0));
    K02 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi0));
    K10 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi1));
    K11 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi1));
    K12 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi1));
    K20 = integral_over_triangle(node_coords, (phi0*u)*Gradient(phi2));
    K21 = integral_over_triangle(node_coords, (phi1*u)*Gradient(phi2));
    K22 = integral_over_triangle(node_coords, (phi2*u)*Gradient(phi2));

    std::vector<std::vector<double>> K = {{K00, K01, K02},
                                          {K10, K11, K12},
                                          {K20, K21, K22}};

    // N matrix assembling

    double N00, N01, N02, N10, N11, N12, N20, N21, N22; 

    N00 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi0));
    N01 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi0));
    N02 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi0));
    N10 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi1));
    N11 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi1));
    N12 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi1));
    N20 = integral_over_triangle(node_coords, ((Gradient(phi0)*u + phi0*Divirgence(u))*u)*Gradient(phi2));
    N21 = integral_over_triangle(node_coords, ((Gradient(phi1)*u + phi1*Divirgence(u))*u)*Gradient(phi2));
    N22 = integral_over_triangle(node_coords, ((Gradient(phi2)*u + phi2*Divirgence(u))*u)*Gradient(phi2));

    std::vector<std::vector<double>> N = {{N00, N01, N02},
                                          {N10, N11, N12},
                                          {N20, N21, N22}}; 

    double alpha = 0.1409714;
    double beta = 0.1160538;
    double gamma = 0.3590284;

    if (step == 1)
    {
        return{M*localmass + (time_step*alpha)*K*localmass - (time_step*time_step*beta)*N*localmass};
    }
    else
    {
        return{M*localmass + time_step*K*localmass_half - (time_step*time_step*gamma)*N*localmass_half};
    }
}