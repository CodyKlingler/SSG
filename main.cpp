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
#include "include/permute.h"
#include "lp_c++.h"



using namespace std;


/* TODO
    -reserve matrix
    -function for probabability of a single vertex
        -may massively speed up all algorithms
        -can use old probability if there is no self reference?
*/

vector<int> hardest_game_switches(0);

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


void get_best_txt(vector<int> &vec){
    ifstream file;
    file.open("hardest_games/best.txt");
    std::string cur_line;

    std::getline(file, cur_line);

    int n = stoi(cur_line);

    std::getline(file, cur_line);
    std::istringstream iss(cur_line);

    for(int i = 0; i<n; i++){
        int x;
        iss >> x;
        vec.push_back(x);
    }

    file.close();
}

void set_best_txt(vector<int> &best_vec){
    ofstream file;
    file.open("hardest_games/best.txt");
    file << best_vec.size() << endl;

    for(auto it = best_vec.begin(); it != best_vec.end(); it++){
        file << *it << '\t';
    } file << endl;

    file.close();
}

void write_hard_game(SSG gg){
    int v = gg.n;
    ofstream file;
    file.open("hardest_games/g" + to_string(v) + ".txt");
    file << gg;
    file.close();
}

int make_min_harder(SSG &gg, int &max_switches){
    //cout << ".." << flush;

    vector<int> vs(gg.min_vertices.begin(), gg.min_vertices.end());
    for(auto it = vs.begin(); it!= vs.end(); it++){
        if(*it < 3)
            vs.erase(it);
    }

    vector<int> val_e(vs.begin(), vs.end());
    val_e.insert(end(val_e), begin(gg.ave_vertices), end(gg.ave_vertices));
    val_e.insert(end(val_e), begin(gg.max_vertices), end(gg.max_vertices));

    vector<int> ts(0);
    for(int i = 1; i<4; i++){
        vs.push_back(i);
    }

    cout << val_e.size() << endl;
    cout << ts.size()<< endl;


    bool increasing;
    do{
        increasing = false;

        vector<int> v_rand(vs.begin(), vs.end());
        std::shuffle(std::begin(v_rand), std::end(v_rand), rng);

        for(auto vit = v_rand.begin(); vit!=v_rand.end(); vit++){
            int v = *vit;
            vertex_type t_og = gg.type[v];
            int e1 = gg.outgoing_edge[v][0];
            int e2 = gg.outgoing_edge[v][1];

            vertex_type type = vertex_type::min;

            vector<int> a_rand(val_e.begin(), val_e.end());
            std::shuffle(std::begin(a_rand), std::end(a_rand), rng);
            for(auto ait = a_rand.begin(); ait!=a_rand.end(); ait++){
                int a = *ait;

                if(a == v)
                    continue;

                vector<int> b_rand(val_e.begin(), val_e.end());
                std::shuffle(std::begin(b_rand), std::end(b_rand), rng);
                for(auto bit = b_rand.begin(); bit!=b_rand.end(); bit++){
                    int b = *bit;

                    if(b == a || b == v)
                        continue;

                    gg.force_vertex(v, type, a, b);
                    
                    int switches = gg.hoffman_karp_n_iterations();

                    if(switches > max_switches){
                        max_switches = switches;
                        //cout << "SWITCHES: " << switches << endl;
                        //cout << gg << endl;
                        increasing = true;

                        t_og = type;
                        e1 = a;
                        e2 = b;
                        if(switches > hardest_game_switches[gg.n]){
                            hardest_game_switches[gg.n] = switches;
                            set_best_txt(hardest_game_switches);
                            write_hard_game(gg);
                        }
                    }
                }
            }
            gg.force_vertex(v, t_og, e1, e2);
        }
    }while(increasing);
    //cout << ",,";
    return max_switches;
}

int make_game_harder_static_n_max(SSG &gg){
    int max_switches = gg.hoffman_karp_n_iterations();

    vector<int> vs = permutation_in_range(0, gg.n-2);
    vector<int> a_rand = permutation_in_range(0, gg.n);

    vector<int> ts = permutation_in_range(1,3);
    vector<int> tm = {3};

    bool increasing;
    do{
        increasing = false;

        std::shuffle(std::begin(vs), std::end(vs), rng);

        for(auto vit = vs.begin(); vit!=vs.end(); vit++){
            int v = *vit;
            vertex_type t_og = gg.type[v];
            int e1 = gg.outgoing_edge[v][0];
            int e2 = gg.outgoing_edge[v][1];
            
            vector<int>* t_vec = (t_og == vertex_type::max) ? &tm : &ts;

            std::shuffle(std::begin(ts), std::end(ts), rng);

            for(auto tit = (*t_vec).begin(); tit!=(*t_vec).end(); tit++){
                int t = *tit;
                vertex_type type = (vertex_type)t;
                std::shuffle(std::begin(a_rand), std::end(a_rand), rng);

                for(auto ait = a_rand.begin(); ait!=a_rand.end(); ait++){
                    int a = *ait;
                    vector<int> b_rand = permutation_in_range(a, gg.n);

                    for(auto bit = b_rand.begin(); bit!=b_rand.end(); bit++){
                        int b = *bit;
                        gg.force_vertex(v, type, a, b);
                        
                        int switches = gg.hoffman_karp_n_iterations();

                        if(switches > max_switches){
                            max_switches = switches;
                            increasing = true;

                            t_og = type;
                            e1 = a;
                            e2 = b;
                            if(switches > hardest_game_switches[gg.n]){
                                hardest_game_switches[gg.n] = switches;
                                set_best_txt(hardest_game_switches);
                                write_hard_game(gg);
                            }
                        }
                    }
                }
            }
            gg.force_vertex(v, t_og, e1, e2);
        }
    }while(increasing);
    return max_switches;
}

int make_game_harder(SSG &gg, int &max_switches){
    //cout << ".." << flush;

    vector<int> vs = permutation_in_range(0, gg.n-2);

    vector<int> a_rand = permutation_in_range(0, gg.n);

    vector<int> ts = permutation_in_range(1,4);

    bool increasing;
    do{
        increasing = false;

        std::shuffle(std::begin(vs), std::end(vs), rng);

        for(auto vit = vs.begin(); vit!=vs.end(); vit++){
            int v = *vit;
            vertex_type t_og = gg.type[v];
            int e1 = gg.outgoing_edge[v][0];
            int e2 = gg.outgoing_edge[v][1];

            std::shuffle(std::begin(ts), std::end(ts), rng);
            for(auto tit = ts.begin(); tit!=ts.end(); tit++){
                int t = *tit;
                vertex_type type = (vertex_type)t;
                std::shuffle(std::begin(a_rand), std::end(a_rand), rng);

                for(auto ait = a_rand.begin(); ait!=a_rand.end(); ait++){
                    int a = *ait;

                    vector<int> b_rand = permutation_in_range(a, gg.n);

                    for(auto bit = b_rand.begin(); bit!=b_rand.end(); bit++){
                        int b = *bit;

                        gg.force_vertex(v, type, a, b);
            
                        int switches = gg.hoffman_karp_n_iterations();

                        if(switches > max_switches){
                            max_switches = switches;
                            //cout << "SWITCHES: " << switches << endl;
                            //cout << gg << endl;
                            increasing = true;

                            t_og = type;
                            e1 = a;
                            e2 = b;
                            if(switches > hardest_game_switches[gg.n]){
                                hardest_game_switches[gg.n] = switches;
                                set_best_txt(hardest_game_switches);
                                write_hard_game(gg);
                            }
                        }
                    }
                }
            }
            gg.force_vertex(v, t_og, e1, e2);
        }
    }while(increasing);
    //cout << ",,";
    return max_switches;
}

int main(int n, char* args[]){
    srand(time(NULL)); std::cout << std::fixed << std::setprecision(4);
    if(hardest_game_switches.size() == 0){
        get_best_txt(hardest_game_switches);
    }


    SSG g = SSG::read_game("folder/unimproved.txt");

    g.hoffman_karp_n_iterations_strat();

    //cout <<= g;
    cout << endl;


    return 0;

    const int n_games = 10000;

    while(true){
        for(int v = 9; v<10; v++){
            cout << v << endl;
            int max = 0;
            for(int n_g = 0; n_g<n_games; n_g++){
                SSG g = SSG::random_game_n_max(3,v);
                
                int prev_hardest = hardest_game_switches[g.n];
                
                int i = g.hoffman_karp_n_iterations();

                if(i > prev_hardest){
                    hardest_game_switches[g.n] = i;
                    set_best_txt(hardest_game_switches);
                    write_hard_game(g);
                    cout << "v: " << v << "  sw: " << i << endl;
                }
            }
        }
    }
    return 0;
}