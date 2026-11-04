#pragma once

#include <array>
#include <cstddef>
#include <functional>

class RootFindingTask
{
public:
    using Function = std::function<double(double)>;

private:
    double left_;
    double right_;
    double tolerance_;
    std::size_t maxIterations_;

    std::array<double, 2> results_{};
    bool hasResults_ = false;

public:
    RootFindingTask(
        double left,
        double right,
        double tolerance,
        std::size_t maxIterations);

    double solveBisection(const Function& f) const;

    double solveNewton(
        const Function& f,
        const Function& df) const;

    double solveSecant(const Function& f) const;

    void solveBoth(
        const Function& f,
        const Function& df);

    bool hasResults() const noexcept;

    const std::array<double, 2>& results() const;
};