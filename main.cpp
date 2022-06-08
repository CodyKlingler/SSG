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
 

 SSG condon_game(){
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

    return a;
 }

std::vector<bool> optimal_condon = {0,1,0,1,1,0,0,0,0,0};

void test_condon_example(){
    SSG game = condon_game();
    auto condon_prob = game.probabilities(optimal_condon);

    vector<vector<bool>> bad_strategies = {{0,0,0,1,0,0,1,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0},
                                            {1,1,1,1,1,1,1,1,1,1},
                                            {0,1,0,1,0,1,0,1,0,1},
                                            {1,0,1,0,1,0,1,0,1,0}};

    vector<vector<bool>> computed_strats(bad_strategies.size());
    vector<vector<double>> computed_probs(bad_strategies.size());

    for(vector<bool>s: bad_strategies){
        vector<bool> hoffman_strat = game.hoffman_karp(s);
        computed_strats.push_back(hoffman_strat);

        auto hoff_prob = game.probabilities(s);
        computed_probs.push_back(hoff_prob);
        
        for(int i = 0; i<condon_prob.size(); i++){
            double delta_prob = condon_prob[i] - hoff_prob[i];
            delta_prob = abs(delta_prob);

            if(delta_prob > .001){
                cout << "TEST_CONDON_EXAMPLE: computed probability does not match optimal." << endl;
                printf("\tvertex: %i\t optimal: %f\t test: %f\t\n", i, condon_prob[i], hoff_prob[i]);
            }
        }
    }

    for(bool b: optimal_condon){
        cout << b << " ";
    }
    //cout << "<-- optimal" << endl;
    for(vector<bool> s: computed_strats){
        for(bool b: s){
            cout << b << " ";
        } 
        cout<<endl << flush;
    }cout << endl;


    cout << flush;
    for(double d: condon_prob){
        printf("%.2f ", d);
        fflush(stdout);
    }
       //cout << "<-- optimal" << endl << flush;
    for(vector<double> s: computed_probs){
        for(double d: s){
            printf("%.2f ", d);
            fflush(stdout);
        } cout<<endl << flush;
    }cout << endl;

    cout << "test of condon example complete." << endl;
}

void test_hoffman_random(int n_tests, int n_vertices){
    for(int i = 0; i<n_tests; i++){
        cout << "testing random game: " << i << "/" << n_tests << endl;
        SSG game = SSG::random_game_loopless(n_vertices);

        const int n_strats_per_game = 10;
        vector<vector<double>> opts(n_strats_per_game);
        vector<vector<bool>> opts_strat(n_strats_per_game);

        for(int k = 0; k<n_strats_per_game; k++){
            cout << "..";
            auto rand_strat = SSG::random_strategy(game.n);
            vector<bool> r_opt = game.hoffman_karp(rand_strat);
            vector<double> opt_prob = game.probabilities(r_opt);
            opts.push_back(opt_prob);
            opts_strat.push_back(r_opt);
        } cout << endl;

        for(int j = 0; j<game.n; j++){
            cout << ",,";
            double amax = 1;
            double amin = -1;
            for(vector<double> cur_opt: opts){
                if(j>=cur_opt.size()){
                    continue;
                }
                amax= (cur_opt[j] > amax)? cur_opt[j] : amax;
                amin = (cur_opt[j] < amin)? cur_opt[j] : amin;
            }

            double delta = amax - amin;
            delta = abs(delta);

            if(delta > .001){
                cout << "test_hoffman_random: something is likely wrong. difference in optimal probabilities greater than .001" << endl;
                cout << delta << endl;
            }
        } cout << endl;

        for(vector<bool> s: opts_strat){
            for(bool b: s){
                cout << b << " ";
            }
            cout << endl;
        }cout << endl;

    }
    cout << "test_hoffman_random complete. "<< endl;
}

int main(){
    srand(time(NULL));

    test_condon_example();

    //test_hoffman_random(5, 15);

    return 0;
}