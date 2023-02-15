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

#include "include/SSG.h"
#include "include/SSG_tests.h"
#include "include/permute.h"

using namespace std;


string cur_folder = "junk"; //folder to write improved games to

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
    file.open(cur_folder + "/best.txt");
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
    file.open(cur_folder + "/best.txt");
    file << best_vec.size() << endl;

    for(auto it = best_vec.begin(); it != best_vec.end(); it++){
        file << *it << '\t';
    } file << endl;

    file.close();
}

void write_hard_game(SSG gg){
    int v = gg.n;
    ofstream file;
    file.open(cur_folder + "/g" + to_string(v) + ".txt");
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

    vector<int> ts = permutation_in_range(2,4);
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

int make_game_harder_static_type_nontrivial(SSG &gg){

    int max_switches = gg.hoffman_karp_n_iterations_inverse();

    vector<int> vs = permutation_in_range(0, gg.n-2);
    vector<int> a_rand = permutation_in_range(0, gg.n);

    vector<int> tm = {3};

    bool increasing;
    do{
        increasing = false;

        std::shuffle(std::begin(vs), std::end(vs), rng);

        for(auto vit = vs.begin(); vit!=vs.end(); vit++){
            int v = *vit;
            int e1 = gg.outgoing_edge[v][0];
            int e2 = gg.outgoing_edge[v][1];


            std::set<int> exclude;
            exclude.insert(v);
            
            for(int i_v: gg.incoming_edges[v]){
                if(gg.type[i_v] == gg.type[v] || abs(gg.type[i_v] - gg.type[v])==2){ //dont let same types (or player) vertices point toward eachother
                    exclude.insert(i_v);
                }
            }

            if(gg.type[v] != vertex_type::ave){ //dont let player vertices point toward sink
                exclude.insert(gg.n-1);
                exclude.insert(gg.n-2);
            }

            std::shuffle(std::begin(a_rand), std::end(a_rand), rng);

            for(auto ait = a_rand.begin(); ait!=a_rand.end(); ait++){
                int a = *ait;

                if(exclude.count(a))
                    continue;


                vector<int> b_rand = permutation_in_range(a, gg.n);

                for(auto bit = b_rand.begin(); bit!=b_rand.end(); bit++){
                    int b = *bit;

                    if(exclude.count(b))
                        continue;
                    if(a == b)
                        continue;

                    gg.force_vertex(v, gg.type[v], a, b);
                    
                    int switches = gg.hoffman_karp_n_iterations_inverse();

                    if(switches > max_switches){
                        max_switches = switches;
                        increasing = true;

                        e1 = a;
                        e2 = b;
                    }
                    if(switches > hardest_game_switches[gg.n]){
                        hardest_game_switches[gg.n] = switches;
                        set_best_txt(hardest_game_switches);
                        write_hard_game(gg);
                    }
                }
            }
            gg.force_vertex(v, gg.type[v], e1, e2);
        }
    }while(increasing);
    return max_switches;
}

int make_game_harder(SSG &gg){

    int max_switches = gg.hoffman_karp_n_iterations();
    vector<int> vs = permutation_in_range(0, gg.n-2);
    vector<int> a_rand = permutation_in_range(0, gg.n);
    vector<int> ts = permutation_in_range(2,5);

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

int make_game_harder(SSG &gg, vector<bool> init_s){
    //cout << ".." << flush;

    int max_switches = gg.hoffman_karp_n_iterations(init_s);

    vector<int> vs = permutation_in_range(0, gg.n-2);

    vector<int> a_rand = permutation_in_range(0, gg.n);

    vector<int> ts = permutation_in_range(2,5);

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
            
                        int switches = gg.hoffman_karp_n_iterations(init_s);

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

int make_game_harder(SSG &gg, int v_start, int v_end){
    //cout << ".." << flush;
    int max_switches = gg.hoffman_karp_n_iterations();

    vector<int> vs = permutation_in_range(v_start, v_end);

    vector<int> a_rand = permutation_in_range(0, gg.n);

    vector<int> ts = permutation_in_range(2,5);

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

int make_game_harder_constant_type(SSG &gg){
    //cout << ".." << flush;
    int max_switches = gg.hoffman_karp_n_iterations();

    vector<int> vs = permutation_in_range(0, gg.n-2);
    vector<int> a_rand = permutation_in_range(0, gg.n);

    bool increasing;
    do{
        increasing = false;

        std::shuffle(std::begin(vs), std::end(vs), rng);

        for(auto vit = vs.begin(); vit!=vs.end(); vit++){
            int v = *vit;
            vertex_type type = gg.type[v];
            int e1 = gg.outgoing_edge[v][0];
            int e2 = gg.outgoing_edge[v][1];

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
            gg.force_vertex(v, type, e1, e2);
        }
    }while(increasing);
    //cout << ",,";
    return max_switches;
}


int main(int n_args, char* args[]){
    srand(time(NULL)); std::cout << std::fixed << std::setprecision(4);
    
    vector<SSG> games;
    for(int i = 0; i<1000000000; i++){
        auto g = SSG::random_game(20);

        auto s_cfb = g.hoffman_karp_min_LP();
        auto s_hk = g.hoffman_karp_max();

        auto p_cfb = g.probabilities(s_cfb);
        auto p_hk = g.probabilities(s_hk);

        if (!SSG::probs_match(p_cfb, p_hk, .01)){
            std::cout << p_cfb << std::endl;
            std::cout << p_hk << std::endl;
        }
    }
    


     // vertex count vs time
    for(int v = 5; v<1000; v++){
        vector<SSG> games;
        for(int i = 0; i<10; i++){
            games.push_back(SSG::random_game(v));
        }
        benchmark_SSG(games);
    }
    return 0;

    // distribution of time
    vector<SSG> games2;
    for(int n_g = 0; n_g < 1000; n_g++){
        int v = 60;
        games2.push_back(SSG::random_game(v));
    }
    benchmark_SSG_distribution2(games2, .5, 30, .5);
    return 0;

}