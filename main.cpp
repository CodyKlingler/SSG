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
#include "lp_c++.h"


using namespace std;
 

/* TODO 
    -investigate Derman's LP for finding optimal min strategy

*/


int main(){
    srand(time(NULL));
    std::cout << std::fixed << std::setprecision(3);

    return main2();

    #define aba

    #ifdef aba
        std::ifstream myfile;
        myfile.open("folder/game_0.txt");
        SSG g = SSG::read_game_file(myfile);
        myfile.close();


        myfile.open("folder/strat_0.txt");
        auto s = SSG::read_strategy_file(myfile);
        myfile.close();

        auto hk = g.hoffman_karp(s);
        auto hki = g.incorrect_hoffman_karp(s);
        auto tr = g.tripathi_hoffman_karp(s);

        auto hkp = g.probabilities(hk);
        auto hkip = g.probabilities(hki);
        auto trp = g.probabilities(tr);

        cout << hk << "\t" << hkp << endl;
        cout << hki << "\t" << hkip << endl;
        cout << tr << "\t" << trp << endl;

        return 0;
    #endif
    
    for(int j = 6; j< 100; j++){
        cout << j << endl;
        for(int i = 1; i< 500; i++){
            if(!test_correctness(1,1,6)){
            cout << "bad!";
            return 0;
            }
        }
    }
    
/*
    for(int i = 10; i<10000; i+=5){
       benchmark_SSG(1,1,i);
    }
*/
   
    return 0;
}