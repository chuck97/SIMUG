#include "vecmath.hpp"

namespace SIMUG
{
    // the sum of two sparse INMOST vectors
    void vec_plus_vec_sparse(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax)
    {
        for (unsigned int i = idmin; i < idmax; ++i)
        {
            res[i] = v1[i] + v2[i];
        }
    }

    // the difference of two sparse INMOST vectors
    void vec_minus_vec_sparse(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax)
    {
        for (unsigned int i = idmin; i < idmax; ++i)
        {
            res[i] = v1[i] - v2[i];
        }
    }

    // vector number multiplication for sparse INMOST vectors
    void vec_mult_num_sparse(INMOST::Sparse::Vector& b, double scale, unsigned int idmin, unsigned int idmax)
    {
        for (unsigned int i = idmin; i < idmax; ++i)
        {
            b[i]*=scale;
        }
    }

    // matrix vector multiplication for sparse INMOST matrix and sparse INMOST vector
    void mat_mult_vec_sparse(INMOST::Sparse::Matrix& M, INMOST::Sparse::Vector& b, INMOST::Sparse::Vector& res)
    {
        INMOST::Solver::OrderInfo info;
        info.PrepareMatrix(M, 0);
        info.PrepareVector(b);
        info.Update(b);
        BARRIER
        M.MatVec(1.0, b, 0.0, res);
        BARRIER
        info.RestoreVector(b);
        info.RestoreMatrix(M);
    }
}