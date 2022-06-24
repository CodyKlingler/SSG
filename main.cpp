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
    -reserve matrix
    -function for probabability of a single vertex
        -may massively speed up all algorithms
        -can use old probability if there is no self reference?
*/


std::vector<std::vector<bool>(SSG::*)(std::vector<bool>)> algorithms = {&SSG::bruteforce, &SSG::ludwig_iterative,  &SSG::hoffman_karp, &SSG::tripathi_hoffman_karp };
std::vector<const char*> algorithm_names = {"ludwig", "bruteforce"};
std::vector<const char*> algorithm_abbrev = {"lw", "bf"};


const int n_strats = 5;
const int n_verts = 3;


double show_files(int max_file){
    std::ifstream myfile;
    myfile.open("folder/game_0.txt");
    SSG g = SSG::read_game_file(myfile);
    myfile.close();

    vector<vector<bool>> strats(0);

    for(int i = 0; i< max_file; i++){
        std::string num = std::to_string(i);
        myfile.open("folder/strat_" +num+".txt");
        auto s = SSG::read_strategy_file(myfile);
        strats.push_back(s);
        myfile.close();
    }
    std::vector<bool> qq;

    for(uint a = 0; a<algorithms.size(); a++){

        auto algo = algorithms[a];
        for(auto s: strats){
            if(s.size()<=0)
                continue;
            auto tr = (g.*algo)(s);
            auto trp = g.exact_probabilities(tr);
            auto trpp = g.probabilities(tr);
            //cout << algorithm_abbrev[a] << ":  " << tr << "\t" << trp << "\t" << trpp << endl;
            qq = s;
        }
    }
    double og_beta = g.beta;
    //printf("%.32f\n", g.beta);
    while(!SSG::probs_match(g.exact_probabilities((g.*algorithms[0])(qq)),g.exact_probabilities((g.*algorithms[1])(qq)),g.tolerance)){
        g.beta *= 1.25;
        if(g.beta > 1){
            cout << "Correct solution not found for any beta value tested." << endl;
            break;
        }
    }
    double min_valid = g.beta;

    while(!SSG::probs_match(g.exact_probabilities((g.*algorithms[0])(qq)),g.exact_probabilities((g.*algorithms[1])(qq)),g.tolerance)){
        g.beta *= 1.25;
        if(g.beta > 1){
            //cout << "." << endl;
            break;
        }
    }

    double max_valid = g.beta;


    printf("invalid >%.8f,\t invalid <%.8f,\t default:%.8f\n", max_valid, min_valid, og_beta);
    return g.beta/og_beta;
}

int find_bad_game(){
    for(int v = n_verts; v<10; v++){
        cout << v << endl;
        for(int j = 0; j< 2; j++){
                int n_games = 1;
                if(!test_correctness(n_games,n_strats,v)){
                    //show_files(n_strats);
                }
        }
    }
    return 0;
}

int main(){
    srand(time(NULL)); std::cout << std::fixed << std::setprecision(3);

    //show_files(n_strats); return 0;

    //gerby(); return 0;


    const int n_games = 100;

/*
    for(int v = 5; v<100; v++){
        cout << v << " vertices" << endl;

        int n_games_failed = 0;
        for(int g = 0; g<n_games; g++){
            if(!test_correctness(1,1,v)){
                show_files(1);
                n_games_failed++;
            }
        }
        cout << "%" << n_games_failed*100/(double)n_games << " failed for " << v << " vertices" << endl;
    }
*/

    for(int v = 5; v<1000; v++){
            benchmark_SSG(n_games,20,v);
    }


    return 0;
}