#include "gtest/gtest.h"
#include "../lusolve.hpp"

TEST(LUSolve, lusolve)
{
    int const kSize = 3;

    std::vector<std::vector<double>> a;
    a.resize(kSize);
    std::fill_n(a.begin(), kSize, std::vector<double>(kSize));
    a[0][0] = 1;
    a[0][1] = 2;
    a[0][2] = 3;
    a[1][0] = 4;
    a[1][1] = 5;
    a[1][2] = 6;
    a[2][0] = 7;
    a[2][1] = 8;
    a[2][2] = 9;

    std::vector<double> b;
    b.resize(kSize);
    b[0] = 10;
    b[1] = 11;
    b[2] = 12;

    std::vector<double> x = LUSolve::lusolve(a, b);

    EXPECT_EQ(kSize, x.size());
    ASSERT_FLOAT_EQ(-76.0 / 3.0, x[0]);
    ASSERT_FLOAT_EQ(125.0 / 3.0, x[1]);
    ASSERT_FLOAT_EQ(-16, x[2]);
}


TEST(LUSolve, error_by_singular_matrix)
{
    int const kSize = 3;

    std::vector<std::vector<double>> a;
    a.resize(kSize);
    std::fill_n(a.begin(), kSize, std::vector<double>(kSize));

    std::vector<double> b;
    b.resize(kSize);
    b[0] = 10;
    b[1] = 11;
    b[2] = 12;

    EXPECT_THROW(LUSolve::lusolve(a, b), std::runtime_error);
}
