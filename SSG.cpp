#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>

#include <fstream>

#include "include/Strategy.h"
#include "include/SSG.h"

//#define SSG_print
#include <Eigen/Dense>
#include <Eigen/Sparse>

const char* vertex_type_names[] = {
    "undefined",
    "ave",
    "min",
    "max",
    "sink_min",
    "sink_max"
};

SSG::SSG(int vertices){
    n = vertices;
    type = std::vector<vertex_type>(n, undefined);

    outgoing_edge = std::vector<std::vector<int>>(n,std::vector<int>(2,0));

    incoming_edges = std::vector<std::vector<int>>(n,std::vector<int>());

    min_vertices = std::vector<int>();
    max_vertices = std::vector<int>();
    ave_vertices = std::vector<int>();
    min_sink_vertex = -1;
    max_sink_vertex = -1;

    std::vector<bool> strat;
    token = -1;

    n_steps_terminate = 1000*n;

}

void SSG::set_vertex_type(int vertex, vertex_type type){
    vertex_type old_type = this->type[vertex];

    switch (old_type)
    {
    case undefined: break;

    default:
        std::cerr << "WARNING: SSG::set_vertex_type can't set vertex " << vertex << ", already set to type " << type << std::endl;
        return;
        break;
    }
    
    switch(type){
        case undefined: return; break;
        case min: min_vertices.emplace_back(vertex); break;
        case max: max_vertices.emplace_back(vertex); break;
        case ave: ave_vertices.emplace_back(vertex); break;
        case sink_min: 
            if( min_sink_vertex > 0){
                std::cerr << "WARNING: SSG::set_vertex_type can't set vertex " << vertex << " to be a min-sink, vertex " << min_sink_vertex << "is already the min-sink" << std::endl;
                return;
            }
            min_sink_vertex = vertex;
            break;
        case sink_max: 
            if( max_sink_vertex > 0){
                    std::cerr << "WARNING: SSG::set_vertex_type can't set vertex " << vertex << " to be a max-sink, vertex " << max_sink_vertex << "is already the max-sink" << std::endl;
                    return;
                } 
            max_sink_vertex = vertex;
            break;
    }

    this->type[vertex] = type;
}

void SSG::start(int starting_vertex, std::vector<bool> combined_strategy){
    token = starting_vertex;
    strat = combined_strategy;
    n_steps = 0;
}

/*
    1: max wins, 0: min wins, -1: no winner yet.
    min is declared the winner after n_steps_terminate steps.
*/
int SSG::step(int steps_to_take){
    #ifdef SSG_print
        const char* cur_vertex_type_name = vertex_type_names[type[token]];
    #endif


    for(int i = 0; i<steps_to_take; i++){
        switch(type[token]){
            case undefined:{
                std::cerr << "WARNING: SSG::step token at vertex " << token << " is undefined" << std::endl;
                return -999;
                break;
            }
            case ave:{
                int next_edge = random()%2;
                int next_token = outgoing_edge[token][next_edge];

                #ifdef SSG_print
                    std::cout << "SSG::step: " << "n:" << n_steps << "   " << token << "(" << cur_vertex_type_name << ")  ->  " << next_token << std::endl;
                #endif

                token = next_token;
                n_steps++;
                break;
            }
            case sink_max:{
                
                #ifdef SSG_print
                    std::cout << "SSG::step: " << "n:" << n_steps << "   " << token << "(" << cur_vertex_type_name << ")  ->  MAX wins by sink" << std::endl;
                #endif
                return 1;
                break;
            }
            case sink_min:{
                #ifdef SSG_print
                    std::cout << "SSG::step: " << "n:" << n_steps << "   " << token << "(" << cur_vertex_type_name << ")  ->  MIN wins by sink" << std::endl;
                #endif
                return 0;
                break;
            }
            default:{
                bool b = strat[token];
                int next_token = outgoing_edge[token][b];

                #ifdef SSG_print
                    std::cout << "SSG::step: " << "n:" << n_steps << "   " << token << "(" << cur_vertex_type_name << ")  ->  " << next_token << std::endl;
                #endif

                token = next_token;
                n_steps++;
                break;
            }

           
        }
        if(n_steps >= n_steps_terminate){
                #ifdef SSG_print
                        std::cout << "SSG::step: " << "n:" << n_steps << "   " << token << "(" << cur_vertex_type_name << ")  ->  MIN declared winner by timeout..  n_steps: " << n_steps <<std::endl;
                #endif
                return 0;
            }
    }
    return -1;
}

int SSG::play(int starting_vertex, std::vector<bool> combined_strategy){
    start(starting_vertex, combined_strategy);
    return step(n_steps_terminate);
}

double SSG::play_n(int starting_vertex, std::vector<bool> combined_strategy, int n_trials){

    long max_wins = 0;

    for(int i = 0; i<n_trials; i++){
        max_wins += play(starting_vertex, combined_strategy);
    }

    return max_wins/(double)n_trials;
}

void SSG::print_graph(){
    
    std::vector<std::vector<int>> vec_list = { min_vertices, max_vertices, ave_vertices};
    std::vector<const char*> vec_names =     { "min_vertices", "max_vertices", "ave_vertices"};

    for(unsigned int i = 0; i< vec_list.size(); i++){

        const std::vector<int> cur_vec = vec_list[i];
        const char* cur_name = vec_names[i];

        std::cout << cur_name << ":  ";
        for(std::vector<int>::const_iterator it = cur_vec.begin(); it !=  cur_vec.end(); it++){
            std:: cout << *it << " " ;
        }
        std::cout << std::endl;
    }
    std::cout << "sink_max_vertex:  " << max_sink_vertex << std::endl;
    std::cout << "sink_min_vertex:  " << min_sink_vertex << std::endl;

    std::cout << "\nincoming edges:" << std::endl;
     for(int i = 0; i < n; i++){

        const std::vector<int> cur_vec = incoming_edges[i];

        std::cout << "(" << i << ") :  ";

        for(std::vector<int>::const_iterator it = cur_vec.begin(); it !=  cur_vec.end(); it++){
            std:: cout << *it << " " ;
        }
        std::cout << std::endl;
    }


    std::cout << "\noutgoing edges: " << std::endl;

    for(int i = 0; i < n; i++){
        vertex_type cur_type = type[i];
        printf(vertex_type_names[cur_type]);
        printf(" (%i) : %i %i\n", i, outgoing_edge[i][0], outgoing_edge[i][1]);
    }
}

std::vector<double> SSG::probabilities(std::vector<bool> combined_strategy){
    using namespace Eigen;

    std::vector<Triplet<double>> coeffs(n*3);

    for(std::vector<int>::iterator it = max_vertices.begin(); it != max_vertices.end(); it++){
        bool b = combined_strategy[*it];
        int p = outgoing_edge[*it][b];

        coeffs.push_back(Triplet<double>(*it,*it,2));
        if(*it != p){
            coeffs.push_back(Triplet<double>(*it,p,-2));
        }
    }

    for(std::vector<int>::iterator it = min_vertices.begin(); it != min_vertices.end(); it++){
        bool b = combined_strategy[*it];
        int p = outgoing_edge[*it][b];

        coeffs.push_back(Triplet<double>(*it,*it,2));
        if(*it != p){
            coeffs.push_back(Triplet<double>(*it,p,-2));
        }
    }

    for(std::vector<int>::iterator it = ave_vertices.begin(); it != ave_vertices.end(); it++){

        int p1 = outgoing_edge[*it][0];
        int p2 = outgoing_edge[*it][1];

        coeffs.push_back(Triplet<double>(*it,*it,2));
        if(*it != p1){
            coeffs.push_back(Triplet<double>(*it,p1,-1));
        }
        if(*it != p2 && p2 != p1){
            coeffs.push_back(Triplet<double>(*it,p2,-1));
        }
        
    }

    /*
    Eigen::SparseMatrix<long double, Eigen::ColMajor> M;
    VectorXld res;
    //filling up M, filling up load
    Eigen::SparseLU< Eigen::SparseMatrix<long double, Eigen::ColMajor> > LUSolver;
    LUSolver.analyzePattern(M);
    LUSolver.factorize(M);
    res = LUSolver.solve(load);
    */

    coeffs.push_back(Triplet<double>(max_sink_vertex,max_sink_vertex,2));

    coeffs.push_back(Triplet<double>(min_sink_vertex,min_sink_vertex,2));
    
    SparseMatrix<double> mat(n,n); 
    mat.setFromTriplets(coeffs.begin(), coeffs.end());  // fill A and b;

    LeastSquaresConjugateGradient<SparseMatrix<double> > solver;

    solver.setTolerance(0);

    mat.makeCompressed();

    solver.compute(mat);

    VectorXd vec = VectorXd::Zero(n);
    vec(max_sink_vertex) = 2;

    //Use the factors to solve the linear system 
    VectorXd x = solver.solve(vec).cwiseAbs(); 
    
    double* probs_temp = x.data();

    int vec_size = x.rows();


    std::cout << mat << std::endl << vec << std::endl;

    std::cout << "err: " << solver.error();
 
    return std::vector<double>(probs_temp, probs_temp + vec_size);
}

std::vector<bool> SSG::hoffman_karp(){
    std::vector<bool> s(n,0);
    return hoffman_karp(s);
}

bool SSG::optimize_min(std::vector<bool> &s){
    auto p = probabilities(s);
    return optimize_min(s,p);
}

bool SSG::optimize_min(std::vector<bool> &s, std::vector<double> probs){
        const double tolerance = .0001;

        bool min_switched_edge;
        bool ever_switched = false;
        do{
            min_switched_edge = false;
            for(auto cur_v: min_vertices){
                
                s[cur_v] = !s[cur_v];
                auto alternate_prob = probabilities(s);
                s[cur_v] = !s[cur_v];

                double p_cur = probs[cur_v];
                double p_other = alternate_prob[cur_v];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > tolerance && p_other < p_cur){
                    s[cur_v] = !s[cur_v]; //switch edge.
                    probs = alternate_prob; //update probability vector

                    min_switched_edge = true;
                    ever_switched = true;
                }
            }
        }while(min_switched_edge);
    return ever_switched;
}


bool SSG::optimize_max(std::vector<bool> &s){
    auto p = probabilities(s);
    return optimize_max(s,p);
}

bool SSG::optimize_max(std::vector<bool> &s, std::vector<double> probs){
        const double tolerance = .0001;

        bool min_switched_edge;
        bool ever_switched = false;
        do{
            min_switched_edge = false;
            for(auto cur_v: min_vertices){
                
                s[cur_v] = !s[cur_v];
                auto alternate_prob = probabilities(s);
                s[cur_v] = !s[cur_v];

                double p_cur = probs[cur_v];
                double p_other = alternate_prob[cur_v];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > tolerance && p_other < p_cur){
                    s[cur_v] = !s[cur_v]; //switch edge.
                    probs = alternate_prob; //update probability vector

                    min_switched_edge = true;
                    ever_switched = true;
                }
            }
        }while(min_switched_edge);
    return ever_switched;
}

std::vector<bool> SSG::hoffman_karp(std::vector<bool> s){
    int iterations_to_convergence = 0;

    
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        bool max_has_switchable_edge = false;
        bool max_vert_was_switched = false;
        auto probs = probabilities(s);
        do{ 
            iterations_to_convergence++;
            vert_switched_any = false;
            for(int cur_v: max_vertices){
                
                int cur_edge = outgoing_edge[cur_v][s[cur_v]];
                int other_edge = outgoing_edge[cur_v][!s[cur_v]];

                double p_cur = probs[cur_edge];
                double p_other = probs[other_edge];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > .001 && p_other > p_cur){
                    max_has_switchable_edge = true;
                    // NOTICE: random()%2 has 1/2 of selecting current edge.
                        vert_switched_any = true;
                        max_vert_was_switched = true;
                        s[cur_v] = !s[cur_v]; //switch edge.
                }
            }
        }while(max_has_switchable_edge && !max_vert_was_switched);

        vert_switched_any = optimize_min(s);
        probs = probabilities(s);

    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}

std::vector<bool> SSG::incorrect_hoffman_karp(){
    std::vector<bool> s(n,0);
    return incorrect_hoffman_karp(s);
}

std::vector<bool> SSG::incorrect_hoffman_karp(std::vector<bool> s){
    
    int iterations_to_convergence = 0;
    int player_loops = 0;


    auto probs = probabilities(s);

    bool vert_switched_any;

    std::vector<int> player_verts;
    do{ // while(vert_switched_any);
        iterations_to_convergence++;
        vert_switched_any = false;
        player_verts = (player_verts == min_vertices)?max_vertices:min_vertices; //switch player vector
        std::vector<std::vector<int>> player_verts_vec = {min_vertices, max_vertices};
        for(auto player_verts: player_verts_vec){

            bool vert_switched_player;
            do{
                vert_switched_player = false;
                player_loops++;
                for(auto it = player_verts.begin(); it!= player_verts.end(); it++){

                   
                    int cur_v = *it;
                    //int cur_edge = outgoing_edge[cur_v][s[cur_v]];
                    //int other_edge = outgoing_edge[cur_v][!s[cur_v]];
                    
                    s[cur_v] = !s[cur_v];
                    auto alternate_prob = probabilities(s);
                    s[cur_v] = !s[cur_v];

                    double p_cur = probs[cur_v];
                    double p_other = alternate_prob[cur_v];

                    double delta = p_other - p_cur;
                    delta = abs(delta);

                    if(delta > .001){
                        if(p_other > p_cur && player_verts == max_vertices){
                            vert_switched_any = true;
                            s[cur_v] = !s[cur_v]; //switch edge.
                            probs = alternate_prob; //update probability vector
                            vert_switched_player = true;
                        }
                        else if (p_other < p_cur && player_verts == min_vertices){
                            vert_switched_any = true;
                            s[cur_v] = !s[cur_v]; //switch edge.
                            probs = alternate_prob; //update probability vector
                            vert_switched_player = true;
                        }
                    }
                }
            }while(vert_switched_player);

        }
    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}

std::vector<bool> SSG::tripathi_hoffman_karp(){
    std::vector<bool> s(n,0);
    return tripathi_hoffman_karp(s);  
}

std::vector<bool> SSG::tripathi_hoffman_karp(std::vector<bool> s){
    int iterations_to_convergence = 0;
    int player_loops = 0;

    
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        bool max_has_switchable_edge = false;
        bool max_vert_was_switched = false;
        auto probs = probabilities(s);
        do{ 
            iterations_to_convergence++;
            vert_switched_any = false;
            for(int cur_v: max_vertices){
                
                int cur_edge = outgoing_edge[cur_v][s[cur_v]];
                int other_edge = outgoing_edge[cur_v][!s[cur_v]];

                double p_cur = probs[cur_edge];
                double p_other = probs[other_edge];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > .001 && p_other > p_cur){
                    max_has_switchable_edge = true;
                    // NOTICE: random()%2 has 1/2 of selecting current edge.
                    if(random()%2){
                        vert_switched_any = true;
                        max_vert_was_switched = true;
                        s[cur_v] = !s[cur_v]; //switch edge.
                    }
                }
            }
        }while(max_has_switchable_edge && !max_vert_was_switched);

        probs = probabilities(s);
        bool min_switched_edge;
        do{
            min_switched_edge = false;
            player_loops++;
            for(auto cur_v: min_vertices){
                
                s[cur_v] = !s[cur_v];
                auto alternate_prob = probabilities(s);
                s[cur_v] = !s[cur_v];

                double p_cur = probs[cur_v];
                double p_other = alternate_prob[cur_v];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > .001 && p_other < p_cur){
                    s[cur_v] = !s[cur_v]; //switch edge.
                    probs = alternate_prob; //update probability vector

                    min_switched_edge = true;
                    vert_switched_any = true;
                }
            }
        }while(min_switched_edge);
    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}

SSG SSG::random_game_loopless(int n){
    SSG game(n);

    if(n<2)
        return game;

    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);

    for(int i = n-3; i>= 0; i--){
        vertex_type type = (vertex_type)(rand()%3+1);
        int j = random()%(n - i) + i;
        int k = random()%(n - i) + i;

        game.set_vertex(i, type, j, k);
    }

    return game;
}

std::vector<bool> SSG::random_strategy(int n){
    std::vector<bool> s(n,0);

    for(unsigned int i = 0; i<s.size(); i++){
        s[i] = random()%2;
    }
    return s;
}

SSG SSG::read_game_file(std::ifstream &file){
    std::string cur_line;

    std::getline(file, cur_line);

    int n = stoi(cur_line);
    SSG game(n);

    for(int i = 0; i<n; i++){
        std::getline(file, cur_line);

        int v, type, e1, e2;
        std::istringstream iss(cur_line);
        iss >> v >> type >> e1 >> e2;
        game.set_vertex(v,(vertex_type)type,e1,e2);
    }
    return game;
}

std::vector<bool> SSG::read_strategy_file(std::ifstream &file){
    std::vector<bool> vec;
    std::string cur_line;

    while(std::getline(file, cur_line, ' ')){
        int n = stoi(cur_line);
        bool b = (bool)n;
        vec.push_back(b);
    }
    return vec;
}

std::ostream& operator<<(std::ostream& os, const SSG &game)
{
    os << game.n << std::endl;
    for(int i = 0; i<game.n; i++){
        os << i << "\t" << game.type[i] << "\t" << game.outgoing_edge[i][0] << "\t" <<game.outgoing_edge[i][1] << std::endl;
    }
    return os;
}


//PRIVATE SETTERS
void SSG::set_edges(int vertex, int e1, int e2){
    if(e1 > e2){
        int e3 = e1;
        e1 = e2;
        e2 = e3;
    }

    //out of range for vertices
    if(e1 >= n || e1 < 0 || e2 >= n || e2 < 0){
        std::cerr << "WARNING SSG::set_edges e1 or e2 out of bounds" << std::endl;
        return;
    }

    if(outgoing_edge[vertex][0] > 0 ){
        std::cerr << "WARNING SSG::set_edges can't set vertex " << vertex << " outgoing edges, they are already set." << std::endl;
        return;
    }


    outgoing_edge[vertex][0] = e1;
    outgoing_edge[vertex][1] = e2;

    incoming_edges[e1].emplace_back(vertex);
    incoming_edges[e2].emplace_back(vertex);
}

void SSG::set_vertex(int vertex, vertex_type type, int e1, int e2){

    (*this).set_vertex_type(vertex, type);
    (*this).set_edges(vertex, e1, e2);
}


// << operator for strategies
std::ostream& operator<<(std::ostream& os, const std::vector<bool> &vec){
    for(auto it = vec.begin(); it!=vec.end(); it++){
        os << *it << " ";
    }
    return os;
}
// << operator for strategies
std::ostream& operator<<(std::ostream& os, const std::vector<double> &vec){
    for(auto it = vec.begin(); it!=vec.end(); it++){
        os << *it << " ";
    }
    return os;
}
