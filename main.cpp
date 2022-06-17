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
    -tripathi
    -function to convert non-stopping game into a stopping game.
*/


std::vector<std::vector<bool>(SSG::*)(std::vector<bool>)> algorithms = {&SSG::tripathi_hoffman_karp, &SSG::hoffman_karp, &SSG::bruteforce };
std::vector<const char*> algorithm_names = {"tripathi", "hoff-karp", "bruteforce"};
std::vector<const char*> algorithm_abbrev = {"tp", "hk", "bf"};


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
    while(!SSG::probs_match(g.exact_probabilities(g.hoffman_karp(qq)),g.exact_probabilities(g.bruteforce(qq)),.0001)){
        g.beta *= 1.1;
        if(g.beta > 1){
            cout << "SOMEFING WONG AS WONG CAN BE" << endl;
            break;
        }
    }
    double min_valid = g.beta;

    while(SSG::probs_match(g.exact_probabilities(g.hoffman_karp(qq)),g.exact_probabilities(g.bruteforce(qq)),.0001)){
        g.beta *= 1.33333;
        if(g.beta > 1){
            //cout << "." << endl;
            break;
        }
    }

    double max_valid = g.beta;


    printf("%.8f\t%.8f\t%.8f\n", max_valid, min_valid, og_beta);
    return g.beta/og_beta;
}

int find_bad_game(){
    for(int v = n_verts; v<100; v++){
        cout << v << endl;
        for(int j = 0; j< 2; j++){
                int n_games = 1;   
                if(!test_correctness(n_games,n_strats,v)){
                    show_files(n_strats);
                }
        }
    }
    return 0;
}

int find_bad_stopping(    int n_verts){
    SSG g = SSG::random_game(n_verts);
    SSG sg = g.stopping_game();

    auto gg = sg.hoffman_karp();
    auto sgg = sg.bruteforce();
    //auto sghk = sg.hoffman_karp()

    //gg = sg.random_strategy();

    auto ggp = sg.exact_probabilities(gg);
    auto sgp = sg.exact_probabilities(sgg);
    bool match = SSG::probs_match(ggp, sgp, .001);

    if(!match){
        for(int i= 0; i < n_verts; i++){
            cout << sgp[i] << " ";
        } cout << endl;
        for(int i= 0; i < n_verts; i++){
            cout << ggp[i]<< " ";
        } cout << endl;
        //cout << ggp << endl;
        //cout << sgp << endl;

        cout << endl;
        cout << g << endl;
        //cout << sg << endl;

        ofstream myfile;
        myfile.open("folder/g.txt");
        myfile << g;
        myfile.close();
        myfile.open("folder/sg.txt");
        myfile << sg;
        myfile.close();
    }

    return match;
}


int main(){
    srand(time(NULL)); std::cout << std::fixed << std::setprecision(3);
    
    //show_files(n_strats); return 0;


    main2(); return 0;

    for(int i = 0; i<10000; i++){
        if(i%10 == 0)
            cout << i << endl;
        if(!find_bad_stopping(6)){
            break;
        }
    } return 0;

    SSG g = SSG::random_game(7);
    SSG sg = g.stopping_game();

    auto gg = g.hoffman_karp();
    auto sgg = sg.hoffman_karp();

    auto ggp = g.exact_probabilities(gg);
    auto sgp = sg.exact_probabilities(sgg);

    cout << ggp << endl;
    cout << sgp << endl;


    return 0;



    ofstream myfile;
    myfile.open("folder/g.txt");
    myfile << g;
    myfile.close();
    myfile.open("folder/sg.txt");
    myfile << sg;
    myfile.close();

    return 0;

    for(int i = 4; i<100; i++){
        benchmark_SSG(1000,1,i);
    }

    find_bad_game(); return 0;


    for(int i= 4; i<10000; i++){
        double minr = 99;
        double ave = 0;
        for(int j = 0; j<10; j++){
            double r = test_stopping_constant(3,3,i);
            minr = std::min(minr, r);
            ave += r;
    //        if(r == init_c)
                return 0;
        }
        printf("v:%i \tmin: %.7f\tave: %.7f\n", i, minr, ave/100.0);
    }
    return 0;

/*
    for(int i = 10; i<10000; i+=5){
       benchmark_SSG(1,1,i);
    }
*/
   
    return 0;
}