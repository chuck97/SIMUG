#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include "defines.hpp"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace SIMUG
{
    // scalar product
    template<typename T>
    T operator* (const std::vector<T>& v1, const std::vector<T>& v2)
    {
        if (v1.size() != v2.size())
            SIMUG_ERR("Can't compute scalar product - different vec sizes!");

        T res = 0;

        for (size_t i = 0; i < v1.size(); i++)
            res+= v1[i]*v2[i];
        return res;
    }

    // matrix mult vector
    template<typename T>
    std::vector<T> operator* (const std::vector<std::vector<T>>& m, const std::vector<T>& v)
    {
        if (m.size() != v.size())
            SIMUG_ERR("Can't multiply matrix on vector with different size!");

        std::vector<T> res(v.size());

        for (size_t i = 0; i < v.size(); ++i)
        {
            res[i] = 0.0;
            for (size_t j = 0; j < v.size(); ++j)
            {
                res[i] += m[i][j]*v[j];
            }
        }
        return res;
    }

    // vector product
    template<typename T>
    std::vector<T> operator% (const std::vector<T>& v1, const std::vector<T>& v2)
    {
        if (v1.size() != v2.size())
        {
            SIMUG_ERR("Can't compute vector product - different vec sizes!");
        }
        else if ((v1.size() != 2) and (v1.size() != 3))
        {
            SIMUG_ERR("Can't compute vector product - available vec sizes: 2, 3!");
        }

        std::vector<T> a = v1;
        std::vector<T> b = v2;

        if (v1.size() == 2)
        {
            a.push_back(0); b.push_back(0);
        }

        return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
    }

    // sum of vectors
    template<typename T>
    std::vector<T> operator+ (const std::vector<T>& v1, const std::vector<T>& v2)
    {
        if (v1.size() != v2.size())
            SIMUG_ERR("Can't compute vectors sum - different vec sizes!");

        std::vector<T> res = v1;

        for (size_t i = 0; i < v1.size(); ++i)
            res[i] = res[i] + v2.at(i);

        return res;
    }

    // difference of vectors
    template<typename T>
    std::vector<T> operator- (const std::vector<T>& v1, const std::vector<T>& v2)
    {
        if (v1.size() != v2.size())
            SIMUG_ERR("Can't compute vectors diff - different vec sizes!");

        std::vector<T> res = v1;

        for (size_t i = 0; i < v1.size(); ++i)
            res[i] = res[i] - v2.at(i);

        return res;
    }

    // vector mult number
    template<typename T>
    std::vector<T> operator* (const std::vector<T>& v, const T& num)
    {
        std::vector<T> res = v;

        for (size_t i = 0; i < v.size(); ++i)
            res[i] *= num;

        return res;
    }

    // number mult vector
    template<typename T>
    std::vector<T> operator* (const T& num, const std::vector<T>& v)
    {
        return v*num;
    }

    // matrix mult number (matrix is stored as vector of raws)
    template<typename T>
    std::vector<std::vector<T>> operator* (const std::vector<std::vector<T>>& m, const T& num)
    {
        std::vector<std::vector<T>> res = m;

        for (size_t i = 0; i < m.size(); ++i)
            res[i] = m[i]*num;

        return res;
    }

    // number mult matrix (matrix is stored as vector of raws)
    template<typename T>
    std::vector<std::vector<T>> operator* (const T& num, const std::vector<std::vector<T>>& m)
    {
        return m*num;
    }

    // L2-norm of vector
    template<typename T>
    T L2_norm_vec(const std::vector<T>& v)
    {
        return std::sqrt(v*v);
    }

    // unit normal for triangle 
    template<typename T>
    std::vector<T> unit_normal(const std::vector<T>& node0, const std::vector<T>& node1, const std::vector<T>& node2)
    {
        if ((node0.size() != node1.size()) or (node0.size() != node2.size()))
        {
            SIMUG_ERR("Can't compute unit normal - different vec coords sizes!");
        }
        else if ((node0.size() != 2) and (node0.size() != 3))
        {
            SIMUG_ERR("Can't compute unit normal - available vec coords sizes: 2, 3!");
        }

        std::vector<T> v01 = node1 - node0;
        std::vector<T> v02 = node2 - node0;

        std::vector<T> normal = v01%v02;
        T normal_size = L2_norm_vec(normal);

        return {normal[0]/normal_size, normal[1]/normal_size, normal[2]/normal_size};
    }

    // angle between two vectors in radians
    template<typename T>
    T angle_vecs(const std::vector<T>& v1, const std::vector<T>& v2)
    {
        if (v1.size() != v2.size())
        {
            SIMUG_ERR("Can't compute angle between two vecs - different vec sizes!");
        }
        else if ((v1.size() != 2) and (v1.size() != 3))
        {
            SIMUG_ERR("Can't compute angle between two vecs - available vec sizes: 2, 3!");
        }

        return acosl(std::min((v1*(1.0/L2_norm_vec(v1)))*(v2*(1.0/L2_norm_vec(v2))), 1.0));
    }

    // square of triangle with given node coords
    template<typename T>
    T trian_square(const std::vector<T>& node0, const std::vector<T>& node1, const std::vector<T>& node2)
    {
        if ((node0.size() != node1.size()) or (node0.size() != node2.size()))
        {
            SIMUG_ERR("Can't compute triangle square - different vec coords sizes!");
        }
        else if ((node0.size() != 2) and (node0.size() != 3))
        {
            SIMUG_ERR("Can't compute triangle square - available vec coords sizes: 2, 3!");
        }

        std::vector<T> v01 = node2 - node0;
        std::vector<T> v02 = node1 - node0;

        return (L2_norm_vec(v01%v02)/2.0);
    }

    // 2x2 or 3x3 matrix determinant (matrix is stored as vector of raws)
    template<typename T>
    T det(const std::vector<std::vector<T>>& m)
    {       
        if (m.size() == 1)
        {
            return m[0][0];
        }
        else if (m.size() == 2)
        {
            T a00 = m[0][0];
            T a01 = m[0][1];
            T a10 = m[1][0];
            T a11 = m[1][1];

            return (a00*a11 - a01*a10);
        }
        else if (m.size() == 3)
        {
            T a00 = m[0][0];
            T a01 = m[0][1];
            T a02 = m[0][2];
            T a10 = m[1][0];
            T a11 = m[1][1];
            T a12 = m[1][2];
            T a20 = m[2][0];
            T a21 = m[2][1];
            T a22 = m[2][2];

            return (a00*a11*a22 + a01*a12*a20 + a10*a21*a02 -
                    a02*a11*a20 - a01*a10*a22 - a12*a21*a00);
        }
        else
        {
            SIMUG_ERR("Can compute determinant only for 1x1, 2x2 or 3x3 matrix!");
        }
    }

    // algebraic addition for 2x2 or 3x3 matrix (matrix is stored as vector of raws)
    template <typename T>
    T alg_add(const std::vector<std::vector<T>>& m, size_t rawnum, size_t colnum)
    {
        if (m.size() > 3)
        {
            SIMUG_ERR("Can compute algebraic addition only for 2x2 or 3x3 matricies");
        }
        else
        {
            std::vector<std::vector<T>> add_matr;
            for (size_t raw = 0; raw < m.size(); ++raw)
            {
                if (raw == rawnum)
                    continue;
                
                std::vector<T> curr_raw;
                
                for (size_t col = 0; col < m.size(); ++col)
                {
                    if (col != colnum)
                        curr_raw.push_back(m[raw][col]);
                }
                add_matr.push_back(curr_raw);
            }
            return (pow(-1, colnum + rawnum)*det(add_matr));
        }
    }
    
    // matrix transposition (matrix is stored as a vector of raws)
    template <typename T>
    std::vector<std::vector<T>> transp(const std::vector<std::vector<T>>& m)
    {
        if (m.size() > 3)
        {
            SIMUG_ERR("Can compute transpose only for 2x2 or 3x3 matricies");
        }
        else
        {
            std::vector<std::vector<T>> res = m;
            
            for (size_t raw = 0; raw < m.size(); ++raw)
            {
                for (size_t col = 0; col < m.size(); ++col)
                {
                    res[raw][col] = m[col][raw];
                }
               
            }
            return res;
        }
    }

    // 2x2 or 3x3 inverse matrix (matrix is stored as vector of raws)
    template <typename T>
    std::vector<std::vector<T>> inv(const std::vector<std::vector<T>>& m)
    {
        if (m.size() > 3)
            SIMUG_ERR("can compute inverse matrix only for 2x2 or 3x3 matricies!");
        
        if (std::abs(det(m)) < REAL_MIN_ABS_VAL)
            SIMUG_ERR("determinant of matrix is almost zero - can't inverse!");

        std::vector<std::vector<T>>  inv_m = m;

        for (size_t raw = 0; raw < m.size(); ++raw)
        {
            for (size_t col = 0; col < m.size(); ++col)
            {
                inv_m[raw][col] = alg_add(transp(m), raw, col);
            }
        }
        return inv_m*(1.0/det(m));
    }

    // output operator for vector
    template <typename T>
    std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
    {
        for (int i = 0; i < v.size(); ++i)
            out << v[i] << " ";
        return out;
    }
    
    // output operator for matrix
    template <typename T>
    std::ostream& operator<< (std::ostream& out, const std::vector<std::vector<T>>& m)
    {
        for (int i = 0; i < m.size(); ++i)
        {
            for (int j = 0; j < m.size(); ++j)
                out << m[i][j] << " ";
            out << std::endl;
        }
        return out;
    }

    // matrix multiplication
    template <typename T>
    std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& lhs, const std::vector<std::vector<T>>& rhs)
    {
        if (lhs.size() != rhs.size())
            SIMUG_ERR("can't multiply square matricies with different sizes!");
            
        std::vector<std::vector<T>> res = lhs;
        for(size_t i = 0; i < rhs.size(); ++i)
            res[i] = transp(rhs)*lhs[i];
        
        return res;
    }

    // solve 2x2 or 3x3 linear system (matrix is stored as a vector of raws)
    template <typename T>
    std::vector<T> solve_linear_system(const std::vector<std::vector<T>>& lhs, const std::vector<T>& rhs)
    {
        if (lhs.size() != rhs.size())
            SIMUG_ERR("can't solve system - LHS matrix and RHS vector should have the same size!");
        
        if ((lhs.size() > 3) or (lhs.size() < 2))
            SIMUG_ERR("can't solve system  - can solve only 2x2 or 3x3 linear systems!");
        
        if (std::abs(det(lhs)) < REAL_MIN_ABS_VAL)
            SIMUG_ERR("can't solve system  - determinant of lhs is to low!");
        
        return (inv(lhs)*rhs);
    }
}