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
    -reserve matrix
    -find function for SSG.c and SSG.tolerance

*/

int find_bad_game(){
    for(int v = 10; v<100; v++){
        cout << v << endl;
        for(int j = 0; j< 500; j++){
                if(!test_correctness(1,10,v)){
                    cout << "bad!";
                    return 0;
                }
        }
    }
    return 0;
}


int main(){
    srand(time(NULL));
    std::cout << std::fixed << std::setprecision(3);
    const double init_c = .0000000001;
    //return main2();

    
    find_bad_game(); return 0;

    for(int i= 4; i<10000; i++){
        double minr = 99;
        double ave = 0;
        for(int j = 0; j<100; j++){
            double r = test_stopping_constant(3,3,i);
            minr = std::min(minr, r);
            ave += r;
            if(r == init_c)
                return 0;
        }
        printf("v:%i \tmin: %.7f\tave: %.7f\n", i, minr, ave/100.0);
    }

    return 0;

    std::cout << endl;
    
    #define aba

    #ifdef aba
        std::ifstream myfile;
        myfile.open("folder/game_0.txt");
        SSG g = SSG::read_game_file(myfile);
        myfile.close();


        int max_file = 9;
        vector<vector<bool>> strats(0);

        for(int i = 0; i<= max_file; i++){
            std::string num = std::to_string(i);
            myfile.open("folder/strat_" +num+".txt");
            auto s = SSG::read_strategy_file(myfile);
            strats.push_back(s);
            myfile.close();
        }

        for(auto s: strats){
            auto hk = g.hoffman_karp(s);
            auto hki = g.incorrect_hoffman_karp(s);
            auto tr = g.tripathi_hoffman_karp(s);

            auto hkp = g.exact_probabilities(hk);
            auto hkip = g.exact_probabilities(hki);
            auto trp = g.exact_probabilities(tr);

            cout << "hk: "<<  hk << "\t" << hkp << endl;
            //cout << hki << "\t" << hkip << endl;
            cout << "tr: " << tr << "\t" << trp << endl;
        }

        while(g.exact_probabilities(g.hoffman_karp(strats[0]))[0]>.5){
            g.c *= 10;
        }

        printf("%.15f\n", g.c);


        return 0;
    #endif
    

    
/*
    for(int i = 10; i<10000; i+=5){
       benchmark_SSG(1,1,i);
    }
*/
   
    return 0;
}