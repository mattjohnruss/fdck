#pragma once

#include <unordered_map>

namespace mjrfd
{
    namespace stencil
    {
        namespace central_1
        {
            const std::unordered_map<int, double> weights = {
                {-1, -0.5}, {1, 0.5}
            };
        }

        namespace forward_1
        {
            const std::unordered_map<int, double> weights = {
                {0, -1.5}, {1, 2.0}, {2, -0.5}
            };
        }

        namespace backward_1
        {
            const std::unordered_map<int, double> weights = {
                {0, 1.5}, {-1, -2.0}, {-2, 0.5}
            };
        }

        namespace central_2
        {
            const std::unordered_map<int, double> weights = {
                {-1, 1.0}, {0, -2.0}, {1, 1.0}
            };
        }

        namespace first_order
        {
            namespace forward_1
            {
                const std::unordered_map<int, double> weights = {
                    {0, -1.0}, {1, 1.0}
                };
            }

            namespace backward_1
            {
                const std::unordered_map<int, double> weights = {
                    {0, 1.0}, {-1, -1.0}
                };
            }
        }
    }
}
