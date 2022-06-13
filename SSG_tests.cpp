#include <chrono>

#include "include/SSG_tests.h"


//#define SSG_TEST_PRINT


const char* opt_strat_condon = "X 1 X X 1 X 0 X X X";
std::vector<bool> optimal_condon = {0,1,0,1,1,0,0,0,0,0};



void printt(std::vector<bool> s){
    for(bool b: s){
        std::cout << b << " ";
    }
}

void printt(std::vector<bool> s, SSG game){
    for(int i = 0; i < s.size(); i++){
        vertex_type t = game.type[i];

        if(t == vertex_type::max || t == vertex_type::min){
            std::cout << s[i] << " ";
        }
        else{
            printf(". ");
        }
    }
}

void printt(std::vector<double> s){
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

void test_condon_example(){
    SSG game = condon_game();
    auto condon_prob = game.probabilities(optimal_condon);
                                         //"X 1 X x 1 X 0 X X X"        
    std::vector<std::vector<bool>> bad_strategies = {{0,0,0,0,0,0,0,0,0,0},
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

    std::vector<std::vector<bool>> computed_strats(bad_strategies.size());

    for(std::vector<bool>s: bad_strategies){
        std::vector<bool> hoffman_strat = game.incorrect_hoffman_karp(s);
        computed_strats.push_back(hoffman_strat);
    }

    #ifdef SSG_TEST_PRINT
        for(std::vector<bool> s: computed_strats){
            if(s.size()){
                printt(s); std::cout << std::endl;
            }
        }
    
        std::cout << opt_strat_condon << std::endl;
        

        for(std::vector<bool> s: computed_strats){
            if(s.size()){
                printt(game.probabilities(s)); std::cout << std::endl;
            }
        }

        printt(condon_prob); std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "test of condon example complete." << std::endl;
    #endif
}

void test_hoffman(int n_tests, int n_strats_per_game, int n_vertices){
    for(int i = 0; i<n_tests; i++){
        #ifdef SSG_TEST_PRINT
            std::cout << "testing random game: " << i+1 << "/" << n_tests << std::endl;
        #endif
        SSG game = SSG::random_game_loopless(n_vertices);

        std::vector<std::vector<bool>> rand_strats(n_strats_per_game);
        std::vector<std::vector<bool>> opts_strat(n_strats_per_game);

        for(int k = 0; k<n_strats_per_game; k++){
            auto rand_strat = SSG::random_strategy(game.n);
            std::vector<bool> r_opt = game.incorrect_hoffman_karp(rand_strat);
            rand_strats.push_back(rand_strat);
            opts_strat.push_back(r_opt);
        } 
        #ifdef SSG_TEST_PRINT
            std::cout << std::endl;
            for(int i = 0; i<rand_strats.size(); i++){
                if(rand_strats[i].size() && opts_strat[i].size()){
                    printt(rand_strats[i], game);
                    std::cout << "\t";
                    printt(opts_strat[i], game);
                    std::cout << "\t";
                    auto p = game.probabilities(opts_strat[i]);
                    printt(p); std::cout << std::endl;
                }
            }std::cout << std::endl;

            for(std::vector<bool> s: opts_strat){
                if(s.size()){
                    
                }
            }std::cout << std::endl;
        #endif
    }
    #ifdef SSG_TEST_PRINT
        std::cout << "test_hoffman complete. "<< std::endl;
    #endif
        
}







void test_randomized_hoffman(int n_tests, int n_strats_per_game, int n_vertices){
    for(int i = 0; i<n_tests; i++){
        #ifdef SSG_TEST_PRINT
            std::cout << "testing random game: " << i+1 << "/" << n_tests << std::endl;
        #endif
        SSG game = SSG::random_game_loopless(n_vertices);

        std::vector<std::vector<bool>> rand_strats(n_strats_per_game);
        std::vector<std::vector<bool>> opts_strat(n_strats_per_game);

        for(int k = 0; k<n_strats_per_game; k++){
            auto rand_strat = SSG::random_strategy(game.n);
            std::vector<bool> r_opt = game.tripathi_hoffman_karp(rand_strat);
            rand_strats.push_back(rand_strat);
            opts_strat.push_back(r_opt);
        } 
        #ifdef SSG_TEST_PRINT
            std::cout << std::endl;
            for(int i = 0; i<rand_strats.size(); i++){
                if(rand_strats[i].size() && opts_strat[i].size()){
                    printt(rand_strats[i], game);
                    std::cout << "\t";
                    printt(opts_strat[i], game);
                    std::cout << "\t";
                    auto p = game.probabilities(opts_strat[i]);
                    printt(p); std::cout << std::endl;
                }
            }std::cout << std::endl;

            for(std::vector<bool> s: opts_strat){
                if(s.size()){
                    
                }
            }std::cout << std::endl;
        #endif
    }
    #ifdef SSG_TEST_PRINT
        std::cout << "test_hoffman complete. "<< std::endl;
    #endif
        
}


std::vector<std::vector<bool>(SSG::*)(std::vector<bool>)> SSG_algorithms= {&SSG::incorrect_hoffman_karp, &SSG::tripathi_hoffman_karp, &SSG::hoffman_karp};

void benchmark_SSG(int n_games, int n_strats_per_game, int n_vertices){
    std::cout << "n: " << n_vertices;

    std::vector<SSG> games(n_games, SSG::random_game_loopless(n_vertices));
    std::vector<std::vector<bool>> random_strategies(n_strats_per_game, SSG::random_strategy(n_vertices));

    for(auto cur_algo: SSG_algorithms){
        auto start = std::chrono::system_clock::now();
        for(SSG cur_game: games){
            for(std::vector<bool> cur_strat: random_strategies){
                std::vector<bool> opt = (cur_game.*cur_algo)(cur_strat);
            }
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << '\t' << elapsed_seconds.count() ;
    }   
    std::cout << std::endl;
}