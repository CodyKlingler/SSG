#include <chrono>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string> 

#include "include/SSG_tests.h"


//#define SSG_TEST_PRINT

std::vector<std::vector<bool>(SSG::*)(std::vector<bool>)> SSG_algorithms = { &SSG::bruteforce, &SSG::ludwig_iterative,  &SSG::hoffman_karp, &SSG::tripathi_hoffman_karp };
std::vector<const char*> SSG_algorithm_names = {"hoff-karp", "ludwig"};

std::vector<std::vector<bool>(SSG::*)(std::vector<bool>)> unused_algorithms = {&SSG::incorrect_hoffman_karp};
std::vector<const char*> unused_names = {"incorrect hoff-karp"};


const char* opt_strat_condon = "X 1 X X 1 X 0 X X X";
std::vector<bool> optimal_condon = {0,1,0,1,1,0,0,0,0,0};



void printt(std::vector<bool> s){
    for(bool b: s){
        std::cout << b << " ";
    }
}

void printt(std::vector<bool> s, SSG game){
    for(uint i = 0; i < s.size(); i++){
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


void benchmark_SSG(int n_games, int n_strats_per_game, int n_vertices){
    std::cout << n_vertices;

    std::vector<SSG> games(n_games, SSG::random_game(n_vertices));
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


bool test_correctness(int n_games, int n_strats_per_game, int n_vertices){
    
    int n_games_written = 0;
    bool strats_written = false;
    bool bad_strat_ever_found = false;

    std::vector<SSG> games(n_games, SSG::random_game(n_vertices));
    std::vector<std::vector<bool>> random_strategies(n_strats_per_game, SSG::random_strategy(n_vertices));
    for(SSG cur_game: games){
            //algo  //strategy       //vertex_p  
        std::vector<std::vector<std::vector<double>>> opt_probs;
        for(auto cur_algo: SSG_algorithms){
            std::vector<std::vector<double>> cur_algo_probs;
            for(std::vector<bool> cur_strat: random_strategies){
                std::vector<bool> opt = (cur_game.*cur_algo)(cur_strat);
                auto cur_p = cur_game.exact_probabilities(opt);
                cur_algo_probs.emplace_back(cur_p);
            }
            opt_probs.emplace_back(cur_algo_probs);
        }

        const double tolerance = games[0].tolerance;

        bool bad_strat_found = false;

        for(uint v = 0; v<opt_probs[0][0].size(); v++){
            double max_prob_v = -999;
            double min_prob_v = 999;
            for(uint a = 0; a<opt_probs.size(); a++){
                double max_prob_s = -999;
                double min_prob_s = 999;
                for(uint s = 0; s<opt_probs[0].size(); s++){
                    double cur_p = opt_probs[a][s][v];

                    max_prob_s = cur_p>max_prob_s? cur_p: max_prob_s;
                    max_prob_v = cur_p>max_prob_v? cur_p: max_prob_v;
                    
                    min_prob_s = cur_p<min_prob_s? cur_p: min_prob_s;
                    min_prob_v = cur_p<min_prob_v? cur_p: min_prob_v;
                }

                double diff_s = abs(max_prob_s - min_prob_s);
                if(diff_s > tolerance){
                   // std::cout << "correctness test for " << SSG_algorithm_names[a] << " produces inconsistent strategies" << std::endl;
                    bad_strat_found = true;
                    bad_strat_ever_found = true;
                    break;
                }
            }


            double diff_v = abs(max_prob_v - min_prob_v);
            if(diff_v > tolerance){
               //std::cout << "correctness test produced inconsistent results between algorithms" << std::endl;
                bad_strat_found = true;
                bad_strat_ever_found = true;
                break;
            }
        }
        // write strategies to file. 
        // write game to file

        if(bad_strat_found){

            std::ofstream myfile;

            if(!strats_written){
                for(int i = 0; i<n_strats_per_game; i++){
                    std::string num_str = std::to_string(i);
                    myfile.open ("folder/strat_" + num_str + ".txt");
                    myfile << random_strategies[i];
                    myfile.close();
                }
            }
            std::string num_str = std::to_string(n_games_written++);
            myfile.open("folder/game_" + num_str + ".txt");
            myfile << cur_game;
            myfile.close();
        }
    }   

    return !bad_strat_ever_found;
}


const double init_c = .0000000001;

double global_cur_beta = 0;

bool probs_match(const std::vector<double> &p1, const std::vector<double> &p2, double tolerance){
    for(uint i = 0; i< std::min(p1.size(), p2.size());i++){
        double delta = std::abs(p1[i] - p2[i]);
        if(delta > tolerance){
            if(global_cur_beta == init_c){
                std::cout << p1 << std::endl << p2 << std::endl;
                }
            return false;
        }
    }
    return true;
}


//create test to compare probabilities for c = 0, vs c = .001 for loopless games

double test_stopping_constant(int n_games, int n_strats_per_game, int n_vertices){
   
    SSG game = SSG::random_game(n_vertices);
    game.beta = 0;

    auto opt_s = game.hoffman_karp();
    auto opt_p = game.exact_probabilities(opt_s);

    for(double cur_beta = init_c; cur_beta < 1; cur_beta+=cur_beta*1.3 ){
        global_cur_beta = cur_beta;
        game.beta = cur_beta;
        auto new_s = game.hoffman_karp();
        auto new_p = game.exact_probabilities(new_s);
        
        if(!probs_match(opt_p, new_p, game.tolerance)){
            if(cur_beta == init_c){

                std::ofstream myfile;


                std::string num_str = std::to_string(42);
                myfile.open("folder/game_" + num_str + ".txt");
                myfile << game;
                myfile.close();

                std::cout << opt_s << std::endl << new_s;
            }
            return cur_beta;
        }
    }
    return 1;
}


int find_bad_stopping();

int find_bad_stopping(int n_verts){
    SSG g = SSG::random_game(n_verts);
    SSG sg = g.stopping_game();

    auto hk_s = sg.hoffman_karp();
    auto bf_s = sg.bruteforce();
    //auto sghk = sg.hoffman_karp()

    //gg = sg.random_strategy();

    auto hk_ep = sg.exact_probabilities(hk_s);
    auto bf_ep = sg.exact_probabilities(bf_s);

    auto hk_p = sg.probabilities(hk_s);
    auto bf_p = sg.probabilities(bf_s);

    bool match = SSG::probs_match(hk_ep, bf_ep, g.tolerance);


    if(!match){
        std::ofstream myfile;
        myfile.open("folder/g.txt");
        myfile << g;
        myfile.close();
        myfile.open("folder/sg.txt");
        myfile << sg;
        myfile.close();

        find_bad_stopping();
    }

    return match;
}

int find_bad_stopping(){

    std::ifstream myfile;

    myfile.open("folder/g.txt");

    SSG g = SSG::read_game_file(myfile);
    int n_verts = g.n;

    SSG sg = g.stopping_game();

    auto hk_s = sg.hoffman_karp();
    auto hkn_s = g.hoffman_karp();
    auto bf_s = sg.bruteforce();
    //auto sghk = sg.hoffman_karp()

    //gg = sg.random_strategy();

    auto hk_ep = sg.exact_probabilities(hk_s);
    auto hkn_ep = g.exact_probabilities(hkn_s);
    auto bf_ep = sg.exact_probabilities(bf_s);

    auto hk_p = sg.probabilities(hk_s);
    auto hkn_p = g.probabilities(hkn_s);
    auto bf_p = sg.probabilities(bf_s);

    bool match = SSG::probs_match(hk_ep, bf_ep, g.tolerance);

 
        std::cout << "bruteforce exact: ";
        for(int i= 0; i < n_verts; i++){
            std::cout << bf_ep[i] << " ";
        } std::cout << std::endl;
        std::cout << "hoff_karp exact:  ";
        for(int i= 0; i < n_verts; i++){
            std::cout << hk_ep[i]<< " ";
        } std::cout << std::endl;
        std::cout << "hk_direct exact:  ";
        for(int i= 0; i < n_verts; i++){
            std::cout << hkn_ep[i]<< " ";
        } std::cout << std::endl;

        std::cout << "bruteforce: \t  ";
        for(int i= 0; i < n_verts; i++){
            std::cout << bf_p[i] << " ";
        } std::cout << std::endl;
        std::cout << "hoff_karp: \t  ";
        for(int i= 0; i < n_verts; i++){
            std::cout << hk_p[i]<< " ";
        } std::cout << std::endl;
        std::cout << "hk_direct: \t  ";
        for(int i= 0; i < n_verts; i++){
            std::cout << hkn_p[i]<< " ";
        } std::cout << std::endl;

        std::cout << bf_s << std::endl;
        std::cout << hk_s << std::endl;

        std::cout << std::endl;
        std::cout << g << std::endl;
        //cout << sg << endl;  

    return match;
}
