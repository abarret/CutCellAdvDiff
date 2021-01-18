#ifndef included_LS_SemiLagrangianAdvIntegrator_inc
#define included_LS_SemiLagrangianAdvIntegrator_inc

#include "LS/SemiLagrangianAdvIntegrator.h"

namespace LS
{
inline double
SemiLagrangianAdvIntegrator::evaluateZSpline(const IBTK::VectorNd x, const int order)
{
    double val = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        val *= ZSpline(x(d), order);
    }
    return val;
}

inline int
SemiLagrangianAdvIntegrator::getSplineWidth(const int order)
{
    return order + 1;
}

inline double
SemiLagrangianAdvIntegrator::ZSpline(double x, const int order)
{
    x = abs(x);
    switch (order)
    {
    case 0:
        if (x < 1.0)
            return 1.0 - x;
        else
            return 0.0;
    case 1:
        if (x < 1.0)
            return 1.0 - 2.5 * x * x + 1.5 * x * x * x;
        else if (x < 2.0)
            return 0.5 * (2.0 - x) * (2.0 - x) * (1.0 - x);
        else
            return 0.0;
    case 2:
        if (x < 1.0)
            return 1.0 - 15.0 / 12.0 * x * x - 35.0 / 12.0 * x * x * x + 63.0 / 12.0 * x * x * x * x -
                   25.0 / 12.0 * x * x * x * x * x;
        else if (x < 2.0)
            return -4.0 + 75.0 / 4.0 * x - 245.0 / 8.0 * x * x + 545.0 / 24.0 * x * x * x - 63.0 / 8.0 * x * x * x * x +
                   25.0 / 24.0 * x * x * x * x * x;
        else if (x < 3.0)
            return 18.0 - 153.0 / 4.0 * x + 255.0 / 8.0 * x * x - 313.0 / 24.0 * x * x * x +
                   21.0 / 8.0 * x * x * x * x - 5.0 / 24.0 * x * x * x * x * x;
        else
            return 0.0;
    default:
        TBOX_ERROR("Unavailable order: " << order << "\n");
        return 0.0;
    }
}

inline double
SemiLagrangianAdvIntegrator::weight(const double r)
{
    return std::exp(-r * r);
}
} // namespace LS
#endif
