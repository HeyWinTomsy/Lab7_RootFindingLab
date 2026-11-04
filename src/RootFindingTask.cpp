#include "RootFindingTask.h"
#include "RootFindingExceptions.h"

#include <cmath>
#include <stdexcept>
#include <utility>

RootFindingTask::RootFindingTask(
    double left,
    double right,
    double tolerance,
    std::size_t maxIterations)
    : left_(left),
      right_(right),
      tolerance_(tolerance),
      maxIterations_(maxIterations)
{
    if (!std::isfinite(left_) ||
        !std::isfinite(right_) ||
        left_ >= right_)
    {
        throw InvalidIntervalException();
    }

    if (!std::isfinite(tolerance_) || tolerance_ <= 0.0)
    {
        throw std::invalid_argument(
            "Tolerance must be finite and positive");
    }

    if (maxIterations_ == 0)
    {
        throw std::invalid_argument(
            "Maximum iterations must be positive");
    }
}

// Метод бісекції
double RootFindingTask::solveBisection(
    const Function& f) const
{
    double a = left_;
    double b = right_;

    double fa = f(a);
    double fb = f(b);

    if (!std::isfinite(fa) || !std::isfinite(fb))
    {
        throw NoConvergenceException();
    }

    if (fa == 0.0)
        return a;

    if (fb == 0.0)
        return b;

    if (std::signbit(fa) == std::signbit(fb))
    {
        throw NoSignChangeException();
    }

    for (std::size_t i = 0; i < maxIterations_; ++i)
    {
        const double mid = a + (b - a) / 2.0;
        const double fm = f(mid);

        if (!std::isfinite(fm))
        {
            throw NoConvergenceException();
        }

        if (fm == 0.0 || (b - a) / 2.0 <= tolerance_)
        {
            return mid;
        }

        if (std::signbit(fa) != std::signbit(fm))
        {
            b = mid;
            fb = fm;
        }
        else
        {
            a = mid;
            fa = fm;
        }
    }

    throw NoConvergenceException();
}

// Метод Ньютона
double RootFindingTask::solveNewton(
    const Function& f,
    const Function& df) const
{
    double x = left_ + (right_ - left_) / 2.0;

    for (std::size_t i = 0; i < maxIterations_; ++i)
    {
        const double fx = f(x);
        const double dfx = df(x);

        if (!std::isfinite(fx) || !std::isfinite(dfx))
        {
            throw NoConvergenceException();
        }

        if (std::abs(fx) <= tolerance_)
        {
            return x;
        }

        if (dfx == 0.0)
        {
            throw ZeroDerivativeException();
        }

        const double next = x - fx / dfx;

        if (!std::isfinite(next))
        {
            throw NoConvergenceException();
        }

        if (std::abs(next - x) <= tolerance_)
        {
            if (std::abs(f(next)) <= tolerance_)
            {
                return next;
            }
        }

        x = next;
    }

    throw NoConvergenceException();
}

// Метод січних
double RootFindingTask::solveSecant(
    const Function& f) const
{
    double x0 = left_;
    double x1 = right_;

    double f0 = f(x0);
    double f1 = f(x1);

    for (std::size_t i = 0; i < maxIterations_; ++i)
    {
        if (!std::isfinite(f0) || !std::isfinite(f1))
        {
            throw NoConvergenceException();
        }

        if (std::abs(f0) <= tolerance_)
            return x0;

        if (std::abs(f1) <= tolerance_)
            return x1;

        const double denominator = f1 - f0;

        if (denominator == 0.0)
        {
            throw ZeroSecantDenominatorException();
        }

        const double x2 = x1 - f1 * (x1 - x0) / denominator;

        if (!std::isfinite(x2))
        {
            throw NoConvergenceException();
        }

        if (std::abs(x2 - x1) <= tolerance_)
        {
            if (std::abs(f(x2)) <= tolerance_)
                return x2;
        }

        x0 = x1;
        f0 = f1;

        x1 = x2;
        f1 = f(x1);
    }

    throw NoConvergenceException();
}

void RootFindingTask::solveBoth(
    const Function& f,
    const Function& df)
{

    const double bisectionResult = solveBisection(f);
    const double newtonResult = solveNewton(f, df);

    results_ = {bisectionResult, newtonResult};
    hasResults_ = true;
}

bool RootFindingTask::hasResults() const noexcept
{
    return hasResults_;
}

const std::array<double, 2>& RootFindingTask::results() const
{
    if (!hasResults_)
    {
        throw std::logic_error(
            "No successful results are available");
    }

    return results_;
}