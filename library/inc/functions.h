#pragma once

#include <cmath>

namespace mjrfd
{
    class DifferentiableFunction
    {
    public:
        virtual double value(const double c) const = 0;
        virtual double deriv(const double c) const = 0;
    };

    class ConstantFunction : public DifferentiableFunction
    {
    public:
        ConstantFunction(const double value) :
            value_(value)
        {
        }

        double value(const double c) const override
        {
            return value_;
        }

        double deriv(const double c) const override
        {
            return 0.0;
        }

    private:
        double value_;
    };

    class HillFunction : public DifferentiableFunction
    {
    public:
        HillFunction(const double a, const double n, const double min = 0.0, const double max = 1.0) :
            a_(a), n_(n), min_(min), max_(max)
        {
        }

        double value(const double c) const override
        {
            return min_ + (max_ - min_)*std::pow(c, n_)/(std::pow(a_, n_) + std::pow(c, n_));
        }

        double deriv(const double c) const override
        {
            return (max_ - min_)*(std::pow(a_, n_)*n_*std::pow(c, n_-1))/std::pow(std::pow(a_, n_) + std::pow(c, n_), 2.0);
        }

    private:
        double a_;
        double n_;
        double min_;
        double max_;
    };
}
