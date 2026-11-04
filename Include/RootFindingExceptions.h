#pragma once

#include <stdexcept>
#include <string>

class RootFindingException : public std::runtime_error
{
public:
    explicit RootFindingException(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

class InvalidIntervalException : public RootFindingException
{
public:
    InvalidIntervalException()
        : RootFindingException(
            "Invalid interval: left must be less than right")
    {
    }
};

class NoSignChangeException : public RootFindingException
{
public:
    NoSignChangeException()
        : RootFindingException(
            "Bisection requires a sign change on the interval")
    {
    }
};

class ZeroDerivativeException : public RootFindingException
{
public:
    ZeroDerivativeException()
        : RootFindingException(
            "Newton method encountered a zero derivative")
    {
    }
};

class ZeroSecantDenominatorException : public RootFindingException
{
public:
    ZeroSecantDenominatorException()
        : RootFindingException(
            "Secant method encountered a zero denominator")
    {
    }
};

class NoConvergenceException : public RootFindingException
{
public:
    NoConvergenceException()
        : RootFindingException(
            "Root finding did not converge")
    {
    }
};