#include <gtest/gtest.h>

#include "RootFindingTask.h"
#include "RootFindingExceptions.h"

#include <cmath>
#include <stdexcept>

// Правильність обчислення трьома методами
TEST(RootFindingTest, FindsSquareRootOfTwo)
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

    EXPECT_NEAR(task.solveBisection(f), std::sqrt(2.0), 1e-7);
    EXPECT_NEAR(task.solveNewton(f, df), std::sqrt(2.0), 1e-7);
    EXPECT_NEAR(task.solveSecant(f), std::sqrt(2.0), 1e-7);
}

// Недопустимий інтервал
TEST(RootFindingTest, RejectsInvalidInterval)
{
    EXPECT_THROW(
        RootFindingTask(2.0, 1.0, 1e-6, 100),
        InvalidIntervalException
    );

    EXPECT_THROW(
        RootFindingTask(1.0, 1.0, 1e-6, 100),
        InvalidIntervalException
    );
}

// Бісекція: немає зміни знака
TEST(RootFindingTest, RejectsIntervalWithoutSignChange)
{
    RootFindingTask task(1.0, 2.0, 1e-6, 100);

    const auto f = [](double x)
    {
        return x * x + 1.0;
    };

    EXPECT_THROW(
        task.solveBisection(f),
        NoSignChangeException
    );
}

// Нульова похідна в методі Ньютона
TEST(RootFindingTest, DetectsZeroDerivative)
{
    RootFindingTask task(-1.0, 1.0, 1e-6, 100);

    auto f = [](double x) {
        return x * x * x + 0.1;
    };

    auto df = [](double x) {
        return 3.0 * x * x;
    };

    EXPECT_THROW(
        task.solveNewton(f, df),
        ZeroDerivativeException
    );
}
// Нульовий знаменник методу січних
TEST(RootFindingTest, DetectsZeroSecantDenominator)
{
    RootFindingTask task(-1.0, 1.0, 1e-6, 100);

    const auto f = [](double x)
    {
        return x * x + 1.0;
    };

    EXPECT_THROW(
        task.solveSecant(f),
        ZeroSecantDenominatorException
    );
}

// Перевищення кількості ітерацій
TEST(RootFindingTest, DetectsNonConvergence)
{
    RootFindingTask task(1.0, 2.0, 1e-12, 1);

    const auto f = [](double x)
    {
        return x * x - 2.0;
    };

    EXPECT_THROW(
        task.solveBisection(f),
        NoConvergenceException
    );
}

// Сильна гарантія: старі результати не втрачаються
TEST(RootFindingTest, PreservesResultsAfterFailedSolveBoth)
{
    RootFindingTask task(1.0, 2.0, 1e-8, 100);

    const auto f1 = [](double x)
    {
        return x * x - 2.0;
    };

    const auto df1 = [](double x)
    {
        return 2.0 * x;
    };

    task.solveBoth(f1, df1);

    const auto oldResults = task.results();

    const auto f2 = [](double x)
    {
        return (x - 1.5) * (x - 1.5) * (x - 1.5) + 0.01;
    };

    const auto df2 = [](double x)
    {
        return 3.0 * (x - 1.5) * (x - 1.5);
    };

    EXPECT_THROW(
        task.solveBoth(f2, df2),
        ZeroDerivativeException
    );

    ASSERT_TRUE(task.hasResults());

    EXPECT_DOUBLE_EQ(task.results()[0], oldResults[0]);
    EXPECT_DOUBLE_EQ(task.results()[1], oldResults[1]);
}

// Не можна отримати результати до успішного обчислення
TEST(RootFindingTest, RejectsAccessBeforeSuccessfulCalculation)
{
    RootFindingTask task(1.0, 2.0, 1e-6, 100);

    EXPECT_FALSE(task.hasResults());

    EXPECT_THROW(
        task.results(),
        std::logic_error
    );
}