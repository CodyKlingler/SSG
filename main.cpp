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
 

void printt(vector<bool> s){
    for(bool b: s){
        cout << b << " ";
    }
}

void printt(vector<bool> s, SSG game){
    for(int i = 0; i < s.size(); i++){
        vertex_type t = game.type[i];

        if(t == vertex_type::max || t == vertex_type::min){
            cout << s[i] << " ";
        }
        else{
            printf(". ");
        }
    }
}

void printt(vector<double> s){
    for(double d: s){
        printf("%.2f ", d);
    }
}

 SSG condon_game(){
    SSG a(10);


    a.set_vertex(0, vertex_type::sink_max, 0, 0);
    a.set_vertex(1, vertex_type::max, 2, 5);
    a.set_vertex(2, vertex_type::ave, 1, 3);
    a.set_vertex(3, vertex_type::max, 3, 4);
    a.set_vertex(4, vertex_type::min, 3, 0);
    a.set_vertex(5, vertex_type::ave, 6, 8);
    a.set_vertex(6, vertex_type::max, 9, 7);
    a.set_vertex(7, vertex_type::ave, 6, 0);
    a.set_vertex(8, vertex_type::ave, 9, 5); 
    a.set_vertex(9, vertex_type::sink_min, 9, 9);

    return a;
 }


const char* opt_strat_condon = "X 1 X X 1 X 0 X X X";
std::vector<bool> optimal_condon = {0,1,0,1,1,0,0,0,0,0};

void test_condon_example(){
    SSG game = condon_game();
    auto condon_prob = game.probabilities(optimal_condon);
                                         //"X 1 X x 1 X 0 X X X"        
    vector<vector<bool>> bad_strategies = {{0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,1,0,0,0},
                                           {0,0,0,0,1,0,0,0,0,0},
                                           {0,0,0,0,1,0,1,0,0,0},
                                           {0,0,0,1,0,0,0,0,0,0},
                                           {0,0,0,1,0,0,1,0,0,0},
                                           {0,0,0,1,1,0,0,0,0,0},
                                           {0,0,0,1,1,0,1,0,0,0},
                                           {0,1,0,0,0,0,0,0,0,0},
                                           {0,1,0,0,0,0,1,0,0,0},
                                           {0,1,0,0,1,0,0,0,0,0},
                                           {0,1,0,0,1,0,1,0,0,0},
                                           {0,1,0,1,0,0,0,0,0,0},
                                           {0,1,0,1,0,0,1,0,0,0},
                                           {0,1,0,1,1,0,0,0,0,0},
                                           {0,1,0,1,1,0,1,0,0,0}};

    vector<vector<bool>> computed_strats(bad_strategies.size());

    for(vector<bool>s: bad_strategies){
        vector<bool> hoffman_strat = game.hoffman_karp(s);
        computed_strats.push_back(hoffman_strat);
    }

    for(vector<bool> s: computed_strats){
        if(s.size()){
            printt(s); cout << endl;
        }
    }
    cout << opt_strat_condon << endl;
    
    for(vector<bool> s: computed_strats){
        if(s.size()){
            printt(game.probabilities(s)); cout << endl;
        }
    }
    printt(condon_prob); cout << endl;

    cout << endl;

    cout << "test of condon example complete." << endl;
}

void test_hoffman_random(int n_tests, int n_vertices){
    for(int i = 0; i<n_tests; i++){
        cout << "testing random game: " << i << "/" << n_tests << endl;
        SSG game = SSG::random_game_loopless(n_vertices);

        const int n_strats_per_game = 10;
        vector<vector<bool>> opts_strat(n_strats_per_game);

        for(int k = 0; k<n_strats_per_game; k++){
            auto rand_strat = SSG::random_strategy(game.n);
            vector<bool> r_opt = game.hoffman_karp(rand_strat);
            opts_strat.push_back(r_opt);
        } cout << endl;

        for(vector<bool> s: opts_strat){
            if(s.size()){
                printt(s, game); cout << endl;
            }
        }cout << endl;

        for(vector<bool> s: opts_strat){
            if(s.size()){
                auto p = game.probabilities(s);
                printt(p); cout << endl;
            }
        }cout << endl;

    }
    cout << "test_hoffman_random complete. "<< endl;
}



int main(){
    srand(time(NULL));

    //test_condon_example();  return 0;

    test_hoffman_random(5, 15); return 0;

    SSG game = condon_game();

    vector<bool> s(10, 1);

    auto new_s = game.hoffman_karp(s);

    printt(new_s); cout << endl;

    auto d = game.probabilities(new_s);

    printt(d); 
    
    
    cout << endl;

    return 0;
}