#include "local_assembling.hpp"

namespace SIMUG
{
    std::vector<std::vector<double>> LocalMassMatrixAssembling(const std::vector<std::vector<double>>& node_coords)
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

    double L_ij_entry(const std::vector<std::vector<double>>& node_coords,
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

    Vec2d LocaReferenceMassMatrixAssembling()
    {
        std::vector<std::vector<double>> reference_trian_coords = 
        {
            {0.0, 0.0},
            {0.0, 1.0},
            {1.0, 0.0}
        };

        std::vector<std::vector<double>> coeffs_reference =
        {
            {-1.0, -1.0, 1.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0}
        };

        ScalarFunction phi0(FuncType::linear, coeffs_reference[0], reference_trian_coords);
        ScalarFunction phi1(FuncType::linear, coeffs_reference[1], reference_trian_coords);
        ScalarFunction phi2(FuncType::linear, coeffs_reference[2], reference_trian_coords);

        std::vector<ScalarFunction> phi_vec = {phi0, phi1, phi2};

        std::vector<double> v1d = {0.0, 0.0, 0.0};

        std::vector<std::vector<double>> M_matrix = {v1d, v1d, v1d};

        // compute M matrix
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                M_matrix[i][j] = integral_over_triangle(reference_trian_coords, phi_vec[j]*phi_vec[i]);
            }
        }
        return M_matrix;
    }

    std::pair<Vec3d, Vec3d> LocaReferenceFirstDerivativesMatrixAssembling()
    {
        std::vector<std::vector<double>> reference_trian_coords = 
        {
            {0.0, 0.0},
            {0.0, 1.0},
            {1.0, 0.0}
        };

        std::vector<std::vector<double>> coeffs_reference =
        {
            {-1.0, -1.0, 1.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0}
        };

        ScalarFunction phi0(FuncType::linear, coeffs_reference[0], reference_trian_coords);
        ScalarFunction phi1(FuncType::linear, coeffs_reference[1], reference_trian_coords);
        ScalarFunction phi2(FuncType::linear, coeffs_reference[2], reference_trian_coords);

        std::vector<ScalarFunction> phi_vec = {phi0, phi1, phi2};

        std::vector<std::vector<ScalarFunction>> dphi_matr = 
        {
            {(d_dx(phi0)), d_dx(phi1), d_dx(phi2)},
            {d_dy(phi0), d_dy(phi1), d_dy(phi2)}
        };

        std::vector<double> v1d = {0.0, 0.0, 0.0};
        std::vector<std::vector<double>> v2d = {v1d, v1d, v1d};

        std::vector<std::vector<std::vector<double>>> DX_matrix = {v2d, v2d, v2d};
        std::vector<std::vector<std::vector<double>>> DY_matrix = {v2d, v2d, v2d};

        // compute DX, DY matricies
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    DX_matrix[i][j][k] = integral_over_triangle(reference_trian_coords, phi_vec[j]*phi_vec[k]*dphi_matr[0][i]);
                    DY_matrix[i][j][k] = integral_over_triangle(reference_trian_coords, phi_vec[j]*phi_vec[k]*dphi_matr[1][i]);
                }
            }
        }
        return {DX_matrix, DY_matrix};
    };

    std::vector<Vec4d> LocaReferenceSecondDerivativesMatrixAssembling()
    {
        std::vector<std::vector<double>> reference_trian_coords = 
        {
            {0.0, 0.0},
            {0.0, 1.0},
            {1.0, 0.0}
        };

        std::vector<std::vector<double>> coeffs_reference =
        {
            {-1.0, -1.0, 1.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0}
        };

        ScalarFunction phi0(FuncType::linear, coeffs_reference[0], reference_trian_coords);
        ScalarFunction phi1(FuncType::linear, coeffs_reference[1], reference_trian_coords);
        ScalarFunction phi2(FuncType::linear, coeffs_reference[2], reference_trian_coords);

        std::vector<ScalarFunction> phi_vec = {phi0, phi1, phi2};

        std::vector<double> v1d = {0.0, 0.0, 0.0};
        std::vector<std::vector<double>> v2d = {v1d, v1d, v1d};
        Vec3d v3d = {v2d, v2d, v2d};
        Vec4d v4d = {v3d, v3d, v3d};

        Vec4d XX_matrix = v4d;
        Vec4d XY_matrix = v4d;
        Vec4d YX_matrix = v4d;
        Vec4d YY_matrix = v4d;

        // compute XX, XY, YX, YY matricies
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    for (int l = 0; l < 3; ++l)
                    {
                        XX_matrix[i][j][k][l] = integral_over_triangle(reference_trian_coords, (d_dx(phi_vec[j])*phi_vec[k] + phi_vec[j]*d_dx(phi_vec[k]))*phi_vec[l]*d_dx(phi_vec[i]));
                        XY_matrix[i][j][k][l] = integral_over_triangle(reference_trian_coords, (d_dx(phi_vec[j])*phi_vec[k] + phi_vec[j]*d_dx(phi_vec[k]))*phi_vec[l]*d_dy(phi_vec[i]));
                        YX_matrix[i][j][k][l] = integral_over_triangle(reference_trian_coords, (d_dy(phi_vec[j])*phi_vec[k] + phi_vec[j]*d_dy(phi_vec[k]))*phi_vec[l]*d_dx(phi_vec[i]));
                        YY_matrix[i][j][k][l] = integral_over_triangle(reference_trian_coords, (d_dy(phi_vec[j])*phi_vec[k] + phi_vec[j]*d_dy(phi_vec[k]))*phi_vec[l]*d_dy(phi_vec[i]));
                    } 
                }
            }
        }
        return {XX_matrix, XY_matrix, YX_matrix, YY_matrix};
    };


    std::vector<double> LocalTG2RhsAssembling(const std::vector<std::vector<double>>& node_coords,
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



        std::vector<std::vector<double>> RhsMatrix = M + time_step*K - (time_step*time_step*0.5)*N;

        return {RhsMatrix*localmass};
    }

    std::vector<double> LocalCG2RhsAssembling(const std::vector<std::vector<double>>& node_coords,
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

    std::vector<double> LocalTTG2RhsAssembling(const std::vector<std::vector<double>>& node_coords,
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


    std::vector<double> LocalTTG3RhsAssembling(const std::vector<std::vector<double>>& node_coords,
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

    std::vector<double> LocalTTG4RhsAssembling(const std::vector<std::vector<double>>& node_coords,
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

    Vec2d FastLocalMassMatrixAssembling(const Vec2d& mass_tensor,
                                        const std::vector<double>& Jacobi_info_vec)
    {
        std::vector<std::vector<double>> Jac_inv = 
        {
            {Jacobi_info_vec[0], Jacobi_info_vec[1]},
            {Jacobi_info_vec[2], Jacobi_info_vec[3]}
        };

        double Jacobian = Jacobi_info_vec[4];

        auto M_matr = mass_tensor;

        std::vector<double> v1d = {0.0, 0.0, 0.0};
        std::vector<std::vector<double>> mass_matr = {v1d, v1d, v1d};

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                mass_matr[i][j] += M_matr[i][j]*Jacobian;
            }
        }

        return mass_matr;
    }

    Vec2d FastLocalFirstDerivMatrixAssembling(const std::pair<Vec3d, Vec3d>& first_deriv_tensors,
                                              const std::vector<double>& u_components,
                                              const std::vector<double>& v_components,
                                              const std::vector<double>& Jacobi_info_vec)
    {
        std::vector<std::vector<double>> Jac_inv = 
        {
            {Jacobi_info_vec[0], Jacobi_info_vec[1]},
            {Jacobi_info_vec[2], Jacobi_info_vec[3]}
        };

        double Jacobian = Jacobi_info_vec[4];

        std::vector<double> u_comp = u_components;
        std::vector<double> v_comp = v_components;

        auto Dx = first_deriv_tensors.first;
        auto Dy = first_deriv_tensors.second;

        std::vector<double> v1d = {0.0, 0.0, 0.0};
        std::vector<std::vector<double>> first_deriv_matr = {v1d, v1d, v1d};

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    first_deriv_matr[i][j] += Dx[i][j][k]*(u_comp[k]*Jac_inv[0][0] + v_comp[k]*Jac_inv[0][1])*Jacobian +
                                              Dy[i][j][k]*(u_comp[k]*Jac_inv[1][0] + v_comp[k]*Jac_inv[1][1])*Jacobian;
                }
            }
        }

        return first_deriv_matr;
    }

    Vec2d FastLocalSecondDerivMatrixAssembling(const std::vector<Vec4d>& second_deriv_tensors,
                                               const std::vector<double>& u_components,
                                               const std::vector<double>& v_components,
                                               const std::vector<double>& Jacobi_info_vec)
    {
        std::vector<std::vector<double>> Jac_inv = 
        {
            {Jacobi_info_vec[0], Jacobi_info_vec[1]},
            {Jacobi_info_vec[2], Jacobi_info_vec[3]}
        };

        double Jacobian = Jacobi_info_vec[4];

        std::vector<double> u_comp = u_components;
        std::vector<double> v_comp = v_components;

        auto DD_matr = second_deriv_tensors;

        std::vector<double> v1d = {0.0, 0.0, 0.0};
        std::vector<std::vector<double>> second_deriv_matr = {v1d, v1d, v1d};

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    for (int l = 0; l < 3; ++l)
                    {
                        second_deriv_matr[i][j] += DD_matr[0][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[0][0]*Jac_inv[0][0] + u_comp[k]*v_comp[l]*Jac_inv[0][0]*Jac_inv[0][1] + v_comp[k]*u_comp[l]*Jac_inv[0][1]*Jac_inv[0][0] + v_comp[k]*v_comp[l]*Jac_inv[0][1]*Jac_inv[0][1])*Jacobian + 
                                                   DD_matr[1][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[0][0]*Jac_inv[1][0] + u_comp[k]*v_comp[l]*Jac_inv[0][0]*Jac_inv[1][1] + v_comp[k]*u_comp[l]*Jac_inv[0][1]*Jac_inv[1][0] + v_comp[k]*v_comp[l]*Jac_inv[0][1]*Jac_inv[1][1])*Jacobian + 
                                                   DD_matr[2][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[1][0]*Jac_inv[0][0] + u_comp[k]*v_comp[l]*Jac_inv[1][0]*Jac_inv[0][1] + v_comp[k]*u_comp[l]*Jac_inv[1][1]*Jac_inv[0][0] + v_comp[k]*v_comp[l]*Jac_inv[1][1]*Jac_inv[0][1])*Jacobian +
                                                   DD_matr[3][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[1][0]*Jac_inv[1][0] + u_comp[k]*v_comp[l]*Jac_inv[1][0]*Jac_inv[1][1] + v_comp[k]*u_comp[l]*Jac_inv[1][1]*Jac_inv[1][0] + v_comp[k]*v_comp[l]*Jac_inv[1][1]*Jac_inv[1][1])*Jacobian; 

                    }
                }
            }
        } 

        return second_deriv_matr;
    }
}