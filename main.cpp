#define EIGEN_MPL2_ONLY

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <random>

#include <Eigen/Dense>

#include "include/Strategy.h" //include Strategy.h before SSG.h
#include "include/SSG.h"

using namespace std;

using namespace Eigen;

int main(){
    srand(time(NULL));

    SSG a(10);

    a.set_vertex(1, vertex_type::max, 2, 5);
    a.set_vertex(2, vertex_type::ave, 1, 3);
    a.set_vertex(3, vertex_type::max, 3, 4);
    a.set_vertex(4, vertex_type::min, 3, 0);
    a.set_vertex(5, vertex_type::ave, 6, 8);
    a.set_vertex(6, vertex_type::max, 9, 7);
    a.set_vertex(7, vertex_type::ave, 6, 0);
    a.set_vertex(8, vertex_type::ave, 9, 5); 
    a.set_vertex(9, vertex_type::sink_min, 9, 9);
    a.set_vertex(0, vertex_type::sink_max, 0, 0);
                //   0 1 2 3 4 5 6 7 8 9
    bool strat[10] = {0,1,0,1,1,0,0,0,0,0};
    bool bad_st[10]= {0,0,0,1,0,0,1,0,0,0};

    Strategy s(a.n, strat);
    Strategy bs(a.n, bad_st);
    
    a.probabilities(bs);

    return 0;
}