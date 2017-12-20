#include <iostream>

#include <Eigen/Dense>

int test()
{
    Eigen::MatrixXd m(2,2);

    m(0,0) = 3.0;
    m(1,0) = 2.5;
    m(0,1) = -1.0;
    m(1,1) = m(1,0) + m(0,1);

    std::cout << m << std::endl;

    return 0;
}

using namespace mjfre;

int main()
{
    const unsigned n_node = 101;

    // set up mesh with 101 nodes between 0 and 1
    Mesh1D mesh(0.0, 1.0, n_node);

    for(i = 1; i < n_node-1; ++i)
    {
        mesh.add_equation()
    }

    Problem problem;
    problem.add_mesh(mesh);
    problem.solve();
    
    //return test();
}
