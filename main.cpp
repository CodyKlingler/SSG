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


vector<int> hardest_game_switches(0);



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

    auto rd = std::random_device{};
    auto rng = std::default_random_engine{rd()};


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


int make_game_harder(SSG &gg, int &max_switches){
    //cout << ".." << flush;
    bool increasing;
    do{
        increasing = false;
        for(int v = 0; v<gg.n-2; v++){
            vertex_type t_og = gg.type[v];
            int e1 = gg.outgoing_edge[v][0];
            int e2 = gg.outgoing_edge[v][1];

            for(int t = 1; t<4; t++){
                vertex_type type = (vertex_type)t;
                for(int a = 0; a<gg.n; a++){
                    for(int b = a; b<gg.n; b++){

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



int make_game_harder_rand(SSG &gg, int &max_switches){
    //cout << ".." << flush;

    vector<int> vs(0);
    for(int i = 0; i<gg.n-2; i++){
        vs.push_back(i);
    }

    vector<int> ts(0);
    for(int i = 1; i<4; i++){
        ts.push_back(i);
    }

    auto rd = std::random_device{};
    auto rng = std::default_random_engine{rd()};


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

            vector<int> t_rand(ts.begin(), ts.end());
            std::shuffle(std::begin(t_rand), std::end(t_rand), rng);
            for(auto tit = t_rand.begin(); tit!=t_rand.end(); tit++){
                int t = *tit;
                vertex_type type = (vertex_type)t;

                vector<int> a_rand(vs.begin(), vs.end());
                std::shuffle(std::begin(a_rand), std::end(a_rand), rng);
                for(auto ait = a_rand.begin(); ait!=a_rand.end(); ait++){
                    int a = *ait;
                    vector<int> b_rand(vs.begin()+a, vs.end());
                    std::shuffle(std::begin(b_rand), std::end(b_rand), rng);
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


void find_hard_game(int v){
    const int n_games = 1000;

    int max = 0;

    for(int n_g = 0; n_g<n_games; n_g++){
        SSG g = SSG::random_game_mod(v);
        
        int i = g.hoffman_karp_n_iterations();

        if(i > max){
            max = i;
            std::ofstream myFile;
            myFile.open("folder/hard_game.txt");
            myFile << g;
            myFile.close();
        }
    }
}


void find_hardest_game(int v){
    const int n_games = 1000;

    int hardest_game_sws = 0;
    for(int k = 0; k<10; k++){

        int max = 0;
        for(int n_g = 0; n_g<n_games; n_g++){
            SSG g = SSG::random_game_mod(v);
            
            int i = g.hoffman_karp_n_iterations();

            if(i > max){
                max = i;

                int q = 0;
                make_game_harder(g, q);

                if(q > hardest_game_sws){
                    hardest_game_sws = q;
                    cout << "SWITCHES: " << q << endl;
                    cout << g << endl;
                }
            }
        }
    }
}

bool next_combination(SSG &g, int v){
    if(v<0)
        return false;

    int e1 = g.outgoing_edge[v][0];
    int e2 = g.outgoing_edge[v][1];

    vertex_type type = g.type[v];

    bool e1_overflow = e1 >= g.n-1;
    bool e2_overflow = e2 >= g.n-1 && e1_overflow;
    bool type_overflow = type == vertex_type::max && e2_overflow;

    e2 = e1_overflow ? e2+1 : e2;
    e2 = e2_overflow ? 0 : e2;

    e1 = e1_overflow ? e2 : e1+1;

    type = e2_overflow? (vertex_type)(((int)type)+1) : type;

    type = type_overflow ? vertex_type::ave : type;

    g.force_vertex(v, type, e1, e2);

    return type_overflow ? next_combination(g, v-1) : true;
}

SSG test_all_games(int n){
    SSG g(n);
    for(int i = 0; i<n-2; i++){
        g.set_vertex(i,vertex_type::ave,0,0);
    }
    g.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);
    g.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);

    int max = 0;
    vector<SSG> max_gs(0, SSG(n));

    do{
        int cur = g.hoffman_karp_n_iterations();

        if(cur > max){
            max = cur;
            max_gs.push_back(g.copy());
        }
    }while(next_combination(g, n-3));

    return max_gs.back();
}


bool next_combination_most(SSG &g, int v){
    if(v<0)
        return false;

    int e1 = g.outgoing_edge[v][0];
    int e2 = g.outgoing_edge[v][1];

    vertex_type type = g.type[v];

    bool e2_overflow = e2 >= g.n-1;
    bool e1_overflow = e1 >= g.n-1 && e2_overflow;
    bool type_overflow = type == vertex_type::max && e1_overflow;

    e1 = e2_overflow ? e1+1 : e1;
    e1 = e1_overflow ? 0 : e1;

    e2 = e2_overflow ? e1 : e2+1;

    type = e1_overflow? (vertex_type)(((int)type)+1) : type;

    type = type_overflow ? vertex_type::ave : type;

    //cout << v << " -> "<< e1 << " " << e2 << endl;
    g.force_vertex(v, type, e1, e2);

    
    if(!type_overflow && type == vertex_type::max){
        if(e1 == e2)
            return next_combination_most(g,v);
    }

    if(!type_overflow && type == vertex_type::ave){
        if(e1 == e2 || e1 == v || e2 == v)
            return next_combination_most(g,v);
    }
    

    return type_overflow ? next_combination_most(g, v-1) : true;
}

SSG test_most_games(int n){
    SSG g(n);
    for(int i = 0; i<n-2; i++){
        g.set_vertex(i,vertex_type::ave,0,0);
    }
    g.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);
    g.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);

    double max = 0;
    vector<SSG> max_gs(0, SSG(n));

    do{
        if(g.max_vertices.size()>0){
            int cur = g.hoffman_karp_n_iterations();

            if(cur/(double)g.max_vertices.size() > max){
                max = cur/(double)g.max_vertices.size();
                max_gs.push_back(g.copy());
                cout << g << cur << endl << endl;
            }
        }
    }while(next_combination_most(g, n-3));

    return max_gs.back();
}


void print_max_iterations(){
      for(int v = 6; v<100; v++){
        int actual_max = 0;
        for(int i = 0; i<3; i++){
            find_hard_game(v);
            std::ifstream myfile;
            myfile.open("folder/hard_game.txt");
            SSG ggg = SSG::read_game_file(myfile);

            int max_iters = 0;
            make_game_harder(ggg, max_iters);
            if(max_iters > actual_max)
                actual_max = max_iters;
        }
        printf("v:%i -> %i \n", v, actual_max);
    }
}

int main(int n, char* args[]){
    srand(time(NULL)); std::cout << std::fixed << std::setprecision(1);
    if(hardest_game_switches.size() == 0){
        get_best_txt(hardest_game_switches);
    }

    for(int i = 4; i<32; i++){
        SSG x = test_most_games(i);
        cout << x;
        cout << i << "-->" << x.hoffman_karp_n_iterations() <<endl<<endl;
    }

    return 0;

    /*

    SSG hg = SSG::read_game("folder/hard_game.txt");

    SSG gg = SSG::read_game("hardest_games/g"+to_string(hg.n)+".txt");

    int ig = gg.hoffman_karp_n_iterations();
    cout << endl << endl;
    int ih = hg.hoffman_karp_n_iterations();
    cout << ig << "  " << ih << endl;

    make_min_harder(hg, ih);
    cout <<= gg; cout << endl << endl;

    cout << ig << "  " << ih << endl;

    return 0;
    
    int n_changes = 1;

    while(true){
        SSG master_gg = SSG::read_game("hardest_games/g10.txt");
        SSG gg = SSG::read_game("hardest_games/g10.txt");
        n_changes++;
        cout << "n_changes: " << n_changes << endl;
        for(int p = 0; p<n_changes*10000; p++){

            vector<int> changed_verts(0);
            for(int ch = 0; ch<n_changes; ch++){
                int v = random()%(gg.n-2);
                changed_verts.push_back(v);

                vertex_type type = (vertex_type)((random()%3)+1);
                int e1 = random()%(gg.n-2);
                int e2 = random()%(gg.n-2);
                gg.force_vertex(v, type, e1, e2);
            }

            int i = gg.hoffman_karp_n_iterations();
            if(i > hardest_game_switches[gg.n]){  
                cout << "wow!";
                make_game_harder(gg, i);

                hardest_game_switches[gg.n] = i;
                set_best_txt(hardest_game_switches);
                write_hard_game(gg);

                cout << "SW: " << i << endl;
                cout << gg << endl;
                n_changes = 1;
                break;
            }

            for(auto it = changed_verts.begin(); it!=changed_verts.end(); it++){
                int old_v = *it;
                vertex_type old_t = master_gg.type[old_v];
                int e1 = master_gg.outgoing_edge[old_v][0];
                int e2 = master_gg.outgoing_edge[old_v][1];
                gg.force_vertex(old_v, old_t, e1, e2);
            }

        }
    }
    */

    const int n_games = 10;
    const int v = 14;

    while(true){
        SSG* max_ssg;
        for(int q = 30; q<32; q++){
            int max = 0;
            cout << q << endl;
            for(int n_g = 0; n_g<n_games; n_g++){
                SSG g = SSG::random_game_one_ave(q);
                
                int i = g.hoffman_karp_n_iterations();

                make_game_harder_rand(g, i);

                if(i >= hardest_game_switches[g.n] && g.max_vertices.size() < g.n-2 && g.min_vertices.size() < g.n-2){
                    hardest_game_switches[g.n] = i;
                    set_best_txt(hardest_game_switches);
                    write_hard_game(g);
                    cout << "v: " << q << "  sw: " << i << endl;
                }
            }
        }
    }
    return 0;
}