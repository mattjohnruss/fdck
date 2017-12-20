#pragma once

#include <vector>

#include <residual.h>

namespace mjrfd
{
    class Solver
    {
    public:
        void solve();


    private:
        std::vector<Residual*> m_residuals;
    };
}
