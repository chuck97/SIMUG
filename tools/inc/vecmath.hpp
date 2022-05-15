#pragma once

#include <vector>
#include <math.h>
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

    //vector product
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
            res[i] += v2[i];

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
            res[i] -= v2[i];

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
}