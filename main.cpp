#define EIGEN_MPL2_ONLY

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <iomanip>
#include <chrono>

#include <Eigen/Dense>

#include "include/Strategy.h" //include Strategy.h before SSG.h
#include "include/SSG.h"
#include "include/SSG_tests.h"


using namespace std;
 

/* TODO 
    -investigate Derman's LP for finding optimal min strategy

*/


int main(){
    srand(time(NULL));
    

    SSG g = condon_game();

    vector<bool> s = {0,1,0,1,1,0,0,0,0,0};
    
    auto p = g.probabilities(s);

    std::cout << std::fixed << std::setprecision(3);
    cout << s << "\t" << p << endl;

    return 0;

    for(int i = 6; i< 1000000; i++){
        cout << i << endl;
        if(!test_correctness(i*i*i*i,i*2,i)){
           cout << "bad!";
           return 0;
        }
    }

    for(int i = 10; i<10000; i+=5){
       benchmark_SSG(i,i,i);
    }

   
    return 0;
}