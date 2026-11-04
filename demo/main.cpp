#include "RootFindingTask.h"
#include "RootFindingExceptions.h"

#include <cmath>
#include <exception>
#include <iostream>

int main()
{
    try
    {
        RootFindingTask task(1.0, 2.0, 1e-8, 100);

        const auto f = [](double x)
        {
            return x * x - 2.0;
        };

        const auto df = [](double x)
        {
            return 2.0 * x;
        };

        task.solveBoth(f, df);

        const auto& result = task.results();

        std::cout << "Bisection root: "
                  << result[0] << '\n';

        std::cout << "Newton root: "
                  << result[1] << '\n';

        std::cout << "Secant root: "
                  << task.solveSecant(f) << '\n';
    }
    catch (const RootFindingException& ex)
    {
        std::cerr << "Root finding error: "
                  << ex.what() << '\n';
        return 1;
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: "
                  << ex.what() << '\n';
        return 1;
    }

    return 0;
}