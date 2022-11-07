#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <tuple>



#include "include/Strategy.h"
#include "include/SSG.h"
#include "include/permute.h"
#include "include/gurobi_c++.h"

//#define SSG_print
#include <Eigen/Dense>
#include <Eigen/Sparse>

GRBEnv env = GRBEnv(true);
bool gurobi_started = false;

// Checks to see if a gurobi instance has been created yet. Starting Gurobi for the first time is slow, so it shoudl only be done once.
void start_gurobi(){
    if(! gurobi_started){
        env.set(GRB_IntParam_OutputFlag,0);
        env.start();
        gurobi_started = true;
    }
}


const char* vertex_type_names[] = {
    "undefined",
    "sink_min",
    "min",
    "ave",
    "max",
    "sink_max"
};

SSG::SSG(int vertices){
    n = vertices;
    type = std::vector<vertex_type>(n, undefined);

    outgoing_edge = std::vector<std::vector<int>>(n,std::vector<int>(2,0));

    incoming_edges = std::vector<std::set<int>>(n,std::set<int>());

    min_vertices = std::vector<int>();
    max_vertices = std::vector<int>();
    ave_vertices = std::vector<int>();
    min_sink_vertex = -1;
    max_sink_vertex = -1;

    std::vector<bool> strat;
    token = -1;

    n_steps_terminate = 1000*n;

    beta = 1/pow(2, n*c);
    beta = .01;
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


void SSG::force_vertex_type(int vertex, vertex_type type){
    vertex_type old_type = this->type[vertex];
    this->type[vertex] = type;

    if(type == old_type){
        return;
    }
    else if(old_type != undefined){
        min_vertices.clear();
        max_vertices.clear();
        ave_vertices.clear();

        for(int i = 0; i<n; i++){

            switch(this->type[i]){
                case undefined: return; break;
                case min: min_vertices.push_back(i); break;
                case max: max_vertices.push_back(i); break;
                case ave: ave_vertices.push_back(i); break;
                default: break;
            }
      }
      //std::cout << min_vertices << std::endl;
      //std::cout << max_vertices << std::endl;
      //std::cout << ave_vertices << std::endl;
      return;
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

        auto cur_vec = vec_list[i];
        const char* cur_name = vec_names[i];

        std::cout << cur_name << ":  ";
        for(auto it = cur_vec.begin(); it !=  cur_vec.end(); it++){
            std:: cout << *it << " " ;
        }
        std::cout << std::endl;
    }
    std::cout << "sink_max_vertex:  " << max_sink_vertex << std::endl;
    std::cout << "sink_min_vertex:  " << min_sink_vertex << std::endl;

    std::cout << "\nincoming edges:" << std::endl;
     for(int i = 0; i < n; i++){

        auto cur_vec = incoming_edges[i];

        std::cout << "(" << i << ") :  ";

        for(auto it = cur_vec.begin(); it !=  cur_vec.end(); it++){
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


SSG SSG::copy(){
    SSG g(n);

    for(int i = 0; i<n; i++){
        g.set_vertex(i, type[i], outgoing_edge[i][0], outgoing_edge[i][1]);
    }
    return g;
}


SSG reduce_helper(SSG g){

}

SSG SSG::reduce(){
    

}



//gv + gv*c*n

SSG SSG::stopping_game(){
    int n_sg_vertices = n + c*n*2*n;
    //std::cout << "sg_n: " << n_sg_vertices << std::endl;
    SSG sg(n_sg_vertices); // n + c*edges... edges = n*2

    for(int gv = 0; gv<n; gv++){
        auto gv_type = type[gv];

        sg.set_vertex(gv, gv_type, n+ c*n*2*gv, n+c*n*gv*2+c*n);
        //std::cout << gv << ":  " << n+ c*n*2*gv << "  "<< n+c*n*gv*2+c*n << std::endl;

        for(int e = 0; e <= 1; e++){
            int chain_start = n+c*n*gv*2 + c*n*e;
            int chain_end = n+c*n*gv*2 + c*n*e + c*n-1;

            for(int i=  chain_start; i<chain_end; i++){
                sg.set_vertex(i, vertex_type::ave, outgoing_edge[gv][e], i+1);
            }
            sg.set_vertex(chain_end, vertex_type::ave, outgoing_edge[gv][e], min_sink_vertex);
        }
    }
    sg.beta = 0;
    return sg;
}




std::vector<double> SSG::probabilities(const std::vector<bool> &combined_strategy){
    using namespace Eigen;

    if(beta <= 0) //TODO. Remove???? some probabilities are incorrect with beta = 0
        beta = 0.0000000001;    

    std::vector<Triplet<double>> coeffs(n*3);

    for(std::vector<int>::iterator it = max_vertices.begin(); it != max_vertices.end(); it++){
        bool b = combined_strategy[*it];
        int p = outgoing_edge[*it][b];

        coeffs.push_back(Triplet<double>(*it,*it,1));
        if(*it != p){
            coeffs.push_back(Triplet<double>(*it,p,-1+beta));
        }
    }

    for(std::vector<int>::iterator it = min_vertices.begin(); it != min_vertices.end(); it++){
        bool b = combined_strategy[*it];
        int p = outgoing_edge[*it][b];

        coeffs.push_back(Triplet<double>(*it,*it,1));
        if(*it != p){
            coeffs.push_back(Triplet<double>(*it,p,-1+beta));
        }
    }

    for(std::vector<int>::iterator it = ave_vertices.begin(); it != ave_vertices.end(); it++){

        int p1 = outgoing_edge[*it][0];
        int p2 = outgoing_edge[*it][1];

        coeffs.push_back(Triplet<double>(*it,*it,1+beta));

        double relation = -.5;
        if(*it == p1 || *it == p2 || p1 == p2){
            relation = -1;
        }

        if(*it != p1){
            coeffs.push_back(Triplet<double>(*it,p1,relation));
        }
        if(*it != p2 && p2 != p1){
            coeffs.push_back(Triplet<double>(*it,p2,relation));
        }
    }

    coeffs.push_back(Triplet<double>(max_sink_vertex,max_sink_vertex,1));
    coeffs.push_back(Triplet<double>(min_sink_vertex,min_sink_vertex,1));
    
    Eigen::SparseMatrix< double, Eigen::ColMajor> mat(n,n); 

    mat.setFromTriplets(coeffs.begin(), coeffs.end());  // fill A and b;

    VectorXd vec = VectorXd::Zero(n);
    vec(max_sink_vertex) = 1;

    Eigen::SparseQR< Eigen::SparseMatrix<double, Eigen::ColMajor> , Eigen::COLAMDOrdering<int>> solver;

    mat.makeCompressed();
    solver.compute(mat);

    VectorXd x = solver.solve(vec).cwiseAbs(); 
    
    double* probs_temp = x.data();

    int vec_size = x.rows();
 
    return std::vector<double>(probs_temp, probs_temp + vec_size);
}

std::vector<double> SSG::exact_probabilities(const std::vector<bool> &combined_strategy){
    double prev_beta = beta;
    beta = 0.000000000001;
    auto p = probabilities(combined_strategy);
    beta = prev_beta;
    return p;
}

bool SSG::optimize_min(std::vector<bool> &s){
    auto p = probabilities(s);
    return optimize_min(s,p);
}

bool SSG::optimize_min(std::vector<bool> &s, std::vector<double> probs){

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

        bool min_switched_edge;
        bool ever_switched = false;
        do{
            min_switched_edge = false;
            for(auto cur_v: max_vertices){
                
                s[cur_v] = !s[cur_v];
                auto alternate_prob = probabilities(s);
                s[cur_v] = !s[cur_v];

                double p_cur = probs[cur_v];
                double p_other = alternate_prob[cur_v];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > tolerance && p_other > p_cur){
                    s[cur_v] = !s[cur_v]; //switch edge.
                    probs = alternate_prob; //update probability vector

                    min_switched_edge = true;
                    ever_switched = true;
                }
            }
        }while(min_switched_edge);
    return ever_switched;
}


bool SSG::increment_min_strat(std::vector<bool> &s){
    for(int m: min_vertices){
        bool prev_state = s[m];

        s[m] = !s[m];
        if(!prev_state){
            return true;
        }
    }
    return false;
}

bool SSG::increment_max_strat(std::vector<bool> &s){
    for(int m: max_vertices){
        bool prev_state = s[m];

        s[m] = !s[m];
        if(!prev_state){
            return true;
        }
    }
    return false;
}


std::vector<bool> SSG::strategy_prefer_own(){
    std::vector<bool> s(n,0);
    for(int i = 0; i<n; i++){
        if(type[i] == vertex_type::max){
            s[i] = 1;
        }
    }
    return s;
}

std::vector<bool> SSG::bruteforce_max(std::vector<bool> s){
    std::vector<bool> best_strat(n,0);
    double best_strat_sum = -9e99;

    for(int m: max_vertices){
        s[m] = false;
    }

    do{
        auto p = exact_probabilities(s);
        double sum_p = 0;
        for(double c_p: p){
            sum_p += c_p;
        }

        if(sum_p > best_strat_sum){
            best_strat_sum = sum_p;
            for(int i = 0; i<n; i++){
                best_strat[i] = s[i];
            }
        }

    }while(increment_max_strat(s));
    return best_strat;
}


std::vector<bool> SSG::bruteforce(){
    std::vector<bool> s(n,0);
    return bruteforce(s);
}

std::vector<bool> SSG::bruteforce(std::vector<bool> s){
    std::vector<bool> best_strat(n,0);
    double best_strat_sum = 9e99;

    for(int m: min_vertices){
        s[m] = false;
    }

    do{
        optimize_max(s);
        //s = bruteforce_max(s);
        auto p = exact_probabilities(s);
        double sum_p = 0;
        for(double c_p: p){
            sum_p += c_p;
        }

        if(sum_p < best_strat_sum){
            best_strat_sum = sum_p;
            for(int i = 0; i<n; i++){
                best_strat[i] = s[i];
            }
        }

    }while(increment_min_strat(s));
    optimize_max(best_strat);
    return best_strat;
}

std::vector<bool> SSG::hoffman_karp(){
    std::vector<bool> s(n,0);
    return hoffman_karp(s);
}
std::vector<bool> SSG::hoffman_karp(std::vector<bool> s){
    if(min_vertices.size() < max_vertices.size()){
        return hoffman_karp_max(s);
    }
    return hoffman_karp_min(s);
}

// Hoffman-karp without LP. Finds optimal strat for min-player and switches max strategy.
// Should be used if |min| < |max|
std::vector<bool> SSG::hoffman_karp_max(){
    std::vector<bool> s(n,0);
    return hoffman_karp_max(s);
}
std::vector<bool> SSG::hoffman_karp_max(std::vector<bool> s){
    int iterations_to_convergence = 0;
    
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        vert_switched_any |= optimize_min(s);
        auto probs = probabilities(s);

        bool max_has_switchable_edge = false;
        bool max_vert_was_switched = false;
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

                if(delta > tolerance && p_other > p_cur){
                    max_has_switchable_edge = true;
                    // NOTICE: random()%2 has 1/2 of selecting current edge.
                    vert_switched_any = true;
                    max_vert_was_switched = true;
                    s[cur_v] = !s[cur_v]; //switch edge.
                }
            }
        }while(max_has_switchable_edge && !max_vert_was_switched);

    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}

std::vector<bool> SSG::hoffman_karp_min(){
    std::vector<bool> s(n,0);
    return hoffman_karp_min(s);
}
std::vector<bool> SSG::hoffman_karp_min(std::vector<bool> s){
    int iterations_to_convergence = 0;
    
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        vert_switched_any |= optimize_max(s);
        auto probs = probabilities(s);

        bool min_has_switchable_edge = false;
        bool min_vert_was_switched = false;
        do{ 
            iterations_to_convergence++;
            vert_switched_any = false;
            for(int cur_v: min_vertices){
                
                int cur_edge = outgoing_edge[cur_v][s[cur_v]];
                int other_edge = outgoing_edge[cur_v][!s[cur_v]];

                double p_cur = probs[cur_edge];
                double p_other = probs[other_edge];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > tolerance && p_other < p_cur){
                    min_has_switchable_edge = true;
                    // NOTICE: random()%2 has 1/2 of selecting current edge.
                    vert_switched_any = true;
                    min_vert_was_switched = true;
                    s[cur_v] = !s[cur_v]; //switch edge.
                }
            }
        }while(min_has_switchable_edge && !min_vert_was_switched);

    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}


std::vector<bool> SSG::hoffman_karp2(){
    std::vector<bool> s(n,0);
    return hoffman_karp2(s);
}

bool SSG::reconstruct_strategy_min(std::vector<bool> &strategy, std::vector<double> p){
    bool ret = false;
    for(int i = 0; i<n; i++){
        if(type[i] == vertex_type::min){
            int j = outgoing_edge[i][0];
            int k = outgoing_edge[i][1];
            bool old_strat = strategy[i];
            strategy[i] = p[j] > p[k]; // selects smaller probability.
            if(strategy[i] != old_strat && p[j]!=p[k]){
                ret = true;
            }
        }
    }
    return ret;
}

bool SSG::reconstruct_strategy_max(std::vector<bool> &strategy, std::vector<double> p){
     bool ret = false;
    for(int i = 0; i<n; i++){
        if(type[i] == vertex_type::max){
            int j = outgoing_edge[i][0];
            int k = outgoing_edge[i][1];
            bool old_strat = strategy[i];
            strategy[i] = p[j] < p[k]; // selects bigger probability.
            if(strategy[i] != old_strat && p[j]!=p[k]){
                ret = true;
            }
        }
    }
    return ret;
}

bool SSG::reconstruct_strategy(std::vector<bool> &strategy, std::vector<double> p){
     bool ret = false;
    for(int i = 0; i<n; i++){
        if(type[i] == vertex_type::max || type[i] == vertex_type::min){
            int j = outgoing_edge[i][0];
            int k = outgoing_edge[i][1];
            bool old_strat = strategy[i];
            strategy[i] = p[j] < p[k]; // selects bigger probability. (assuming `i` is a max vertex)
            strategy[i] = type[i] == vertex_type::max ? strategy[i] : !strategy[i]; // flips strategy if is min vertex instead
            if(strategy[i] != old_strat && p[j]!=p[k]){
                ret = true;
            }
        }
    }
    return ret;
}


bool SSG::switch_max(std::vector<bool> &strategy, std::vector<double> p){
    bool any_switched = false;
    for(int i: max_vertices){
        int j = outgoing_edge[i][0];
        int k = outgoing_edge[i][1];
        bool old_strat = strategy[i];
        strategy[i] = p[j] < p[k]; // selects bigger probability.
        if(strategy[i] != old_strat && p[j]!=p[k]){
            any_switched = true;
        }
    }
    return any_switched;
}

bool SSG::switch_min(std::vector<bool> &strategy, std::vector<double> p){
    bool any_switched = false;
    for(int i: min_vertices){
        int j = outgoing_edge[i][0];
        int k = outgoing_edge[i][1];
        bool old_strat = strategy[i];
        strategy[i] = p[j] > p[k]; // selects smaller probability.
        if(strategy[i] != old_strat && p[j]!=p[k]){
            any_switched = true;
        }
    }
    return any_switched;
}

bool SSG::switch_max_tripathi(std::vector<bool> &strategy, std::vector<double> p){
    bool any_switched = false;
    bool any_switchable = false;
    do{
        for(int i: max_vertices){
            int j = outgoing_edge[i][0];
            int k = outgoing_edge[i][1];
            bool old_strat = strategy[i];
            bool optimal_strat = p[j] < p[k]; // selects bigger probability.                    
            if(optimal_strat != old_strat && p[j]!=p[k]){
                any_switchable = true;
                if(random()%2){
                    strategy[i] = optimal_strat;
                    any_switched = true;
                }
            }
        }
    }while(any_switchable && !any_switched);
    return any_switched;
}

std::vector<double> SSG::optimize_min_LP(std::vector<bool> &strategy){
  // returns prob
  std::vector<double> prob(0);

  try {
    start_gurobi();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar vert[n];
    for(int i = 0; i<n; i++){
      vert[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    }

    // Add constraints
    for(int i = 0; i<n; i++){
      int j = outgoing_edge[i][0];
      int k = outgoing_edge[i][1];
      int selected = outgoing_edge[i][strategy[i]];

      double prob_move_next = (1-beta);
        switch(type[i]){
            case vertex_type::max:{
                model.addConstr(vert[i],GRB_EQUAL,vert[selected]*prob_move_next);
                break;
            }
            case vertex_type::min:{
                model.addConstr(vert[i],GRB_LESS_EQUAL,vert[j]*prob_move_next );
                model.addConstr(vert[i],GRB_LESS_EQUAL,vert[k]*prob_move_next );
                break;
            }
            case vertex_type::ave:{
                model.addConstr(vert[i],GRB_EQUAL,.5 *prob_move_next * vert[k] + .5 * prob_move_next  * vert[j]);
                break;
            }
            case vertex_type::sink_max:{
                model.addConstr(vert[i],GRB_EQUAL,1);
                break;
            }
            case vertex_type::sink_min:{
                model.addConstr(vert[i],GRB_EQUAL,0);
                break;
            }
        }
    }

    // find coefficient for each vertex in the objective function.
    std::vector<double> coeffs(n, 1.0);

    // set objective function
    GRBLinExpr objective_fn;
    objective_fn.addTerms(coeffs.data(),vert,n);
    model.setObjective(objective_fn, GRB_MAXIMIZE);

    // solve model
    model.optimize();

    // build probability vector
    for(int i = 0; i<n; i++){
      prob.push_back(vert[i].get(GRB_DoubleAttr_X));
    }
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    std::cout << "Exception during optimization in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return prob;
}

std::vector<double> SSG::optimize_max_LP(std::vector<bool> &strategy){
  // returns prob
  std::vector<double> prob(0);

  try {
    start_gurobi();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar vert[n];
    for(int i = 0; i<n; i++){
      vert[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    }

    // Add constraints
    for(int i = 0; i<n; i++){
      int j = outgoing_edge[i][0];
      int k = outgoing_edge[i][1];
      int selected = outgoing_edge[i][strategy[i]];

      double prob_move_next = (1-beta);
        switch(type[i]){
            case vertex_type::max:{
                model.addConstr(vert[i],GRB_GREATER_EQUAL,vert[j]*prob_move_next );
                model.addConstr(vert[i],GRB_GREATER_EQUAL,vert[k]*prob_move_next );
                break;
            }
            case vertex_type::min:{
                model.addConstr(vert[i],GRB_EQUAL,vert[selected]*prob_move_next);
                break;
            }
            case vertex_type::ave:{
                model.addConstr(vert[i],GRB_EQUAL,.5 *prob_move_next * vert[k] + .5 * prob_move_next  * vert[j]);
                break;
            }
            case vertex_type::sink_max:{
                model.addConstr(vert[i],GRB_EQUAL,1);
                break;
            }
            case vertex_type::sink_min:{
                model.addConstr(vert[i],GRB_EQUAL,0);
                break;
            }
        }
    }

    // find coefficient for each vertex in the objective function.
    std::vector<double> coeffs(n, 1.0);

    // set objective function
    GRBLinExpr objective_fn;
    objective_fn.addTerms(coeffs.data(),vert,n);
    model.setObjective(objective_fn, GRB_MINIMIZE);

    // solve model
    model.optimize();

    // build probability vector
    for(int i = 0; i<n; i++){
      prob.push_back(vert[i].get(GRB_DoubleAttr_X));
    }
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    std::cout << "Exception during optimization in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return prob;
}

//computes optimal strategy for max given a probability vector p.
void SSG::make_max_stable(std::vector<double> &p){
    bool any_changed;
    do{
        
    any_changed = false;

    bool any_changed2;
    do{
        any_changed2 = false;
        for(int x: max_vertices){
            int j = outgoing_edge[x][0];
            int k = outgoing_edge[x][1];
            int max_v = p[j] > p[j] ? j : k;
            if(abs(p[max_v] - p[x]) >= .0000000001){
                p[x] = p[max_v];
                any_changed = true;
                any_changed2 = true;
            }
        }
    }while(any_changed2);

    if(any_changed){
        for(int a: ave_vertices){
            int j = outgoing_edge[a][0];
            int k = outgoing_edge[a][1];
            p[a] = (p[j]+p[k])/2.0;
        }
    }
    }while(any_changed);
}

//computes optimal strategy for max given a probability vector v
std::vector<double> SSG::make_max_stable_LP(std::vector<double> &v){
 // returns prob
  std::vector<double> prob(0);

  try {
    start_gurobi();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar vert[n]; // equivalent to x(i) in condon's paper "On Algorithms for Simple Stochastic Games"
    for(int i = 0; i<n; i++){
      vert[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    }

    // Add constraints
    for(int i = 0; i<n; i++){
        int j = outgoing_edge[i][0];
        int k = outgoing_edge[i][1];

        double prob_move_next = (1-beta); // use (1-beta) to make this a 1/2^m stopping-game

        switch(type[i]){
            case vertex_type::max:{
                    model.addConstr(vert[i],GRB_GREATER_EQUAL,vert[j]*prob_move_next);
                    model.addConstr(vert[i],GRB_GREATER_EQUAL,vert[k]*prob_move_next);
                break;
            }
            case vertex_type::min:{
                if(abs(v[i] - v[j]) < abs(v[i] - v[k])){
                    model.addConstr(vert[i],GRB_EQUAL,vert[j]*prob_move_next);
                }
                else{
                    model.addConstr(vert[i],GRB_EQUAL,vert[k]*prob_move_next);
                }
                
                break;
            }
            case vertex_type::ave:{
                model.addConstr(vert[i],GRB_EQUAL,.5 *prob_move_next * vert[k] + .5 * prob_move_next  * vert[j]);
                break;
            }
            case vertex_type::sink_max:{
                model.addConstr(vert[i],GRB_EQUAL,1);
                break;
            }
            case vertex_type::sink_min:{
                model.addConstr(vert[i],GRB_EQUAL,0);
                break;
            }
        }
    }

    // find coefficient for each vertex in the objective function.
    std::vector<double> coeffs(n, 1.0);

    // set objective function
    GRBLinExpr objective_fn;
    objective_fn.addTerms(coeffs.data(),vert,n);
    model.setObjective(objective_fn, GRB_MINIMIZE);

    // solve model
    model.optimize();

    // build probability vector
    for(int i = 0; i<n; i++){
      prob.push_back(vert[i].get(GRB_DoubleAttr_X));
    }
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    std::cout << "Exception during optimization in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return prob;
}

std::vector<bool> SSG::converge_from_below(){
     //initial `v` vector is all 0's (hence converge from below)

    std::vector<double> current_v(n,0); //TODO: may want to change this to use LP.
    current_v = make_max_stable_LP(current_v);
    std::vector<double> prev_v;
    
    do{
        prev_v = current_v;
        current_v = converge_from_below_LP(prev_v);
        current_v = make_max_stable_LP(current_v);
    }while(!probs_match(prev_v, current_v,tolerance));
    
    std::vector<bool> s(n,0);
    reconstruct_strategy(s,current_v);
    return s;
}

// Input: current probability vector `v`
std::vector<double> SSG::converge_from_below_LP(std::vector<double> v){

    // returns prob
  std::vector<double> prob(0);

  try {
    start_gurobi();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar vert[n]; // equivalent to x(i) in condon's paper "On Algorithms for Simple Stochastic Games"
    for(int i = 0; i<n; i++){
      vert[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    }

    // Add constraints
    for(int i = 0; i<n; i++){
        int j = outgoing_edge[i][0];
        int k = outgoing_edge[i][1];

        double prob_move_next = (1-beta); // use (1-beta) to make this a 1/2^m stopping-game

        switch(type[i]){
            case vertex_type::max:{
                if( abs(v[k] - v[j]) <= tolerance){
                    model.addConstr(vert[i],GRB_EQUAL,.5 *prob_move_next * vert[k] + .5 * prob_move_next  * vert[j]);
                }
                else{
                    int max_v = v[k] > v[j] ? k : j;
                    model.addConstr(vert[i],GRB_EQUAL,vert[max_v]*prob_move_next);
                }
                break;
            }
            case vertex_type::min:{
                model.addConstr(vert[i],GRB_LESS_EQUAL,vert[j]*prob_move_next );
                model.addConstr(vert[i],GRB_LESS_EQUAL,vert[k]*prob_move_next );
                break;
            }
            case vertex_type::ave:{
                model.addConstr(vert[i],GRB_EQUAL,.5 *prob_move_next * vert[k] + .5 * prob_move_next  * vert[j]);
                break;
            }
            case vertex_type::sink_max:{
                model.addConstr(vert[i],GRB_EQUAL,1);
                break;
            }
            case vertex_type::sink_min:{
                model.addConstr(vert[i],GRB_EQUAL,0);
                break;
            }
        }
    }

    // find coefficient for each vertex in the objective function.
    std::vector<double> coeffs(n, 1.0);

    // set objective function
    GRBLinExpr objective_fn;
    objective_fn.addTerms(coeffs.data(),vert,n);
    model.setObjective(objective_fn, GRB_MAXIMIZE);

    // solve model
    model.optimize();

    // build probability vector
    for(int i = 0; i<n; i++){
      prob.push_back(vert[i].get(GRB_DoubleAttr_X));
    }
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    std::cout << "Exception during optimization in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return prob;
}

// non-convex quadratic program to solve the SSG. Often infeasible for n >= 10
std::vector<double> SSG::QP(){
    // returns prob
  std::vector<double> prob(0);

  try {
    start_gurobi();

    // Create an empty model
    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_NonConvex, 2);

    // Create Objective Function
    GRBQuadExpr objective_fn;

    // Create variables
    GRBVar vert[n];
    for(int i = 0; i<n; i++){
      vert[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    }

    // Add constraints
    for(int i = 0; i<n; i++){
      int j = outgoing_edge[i][0];
      int k = outgoing_edge[i][1];

      if(type[i] == vertex_type::max || type[i] == vertex_type::min){
        objective_fn += (vert[i] - vert[j])*(vert[i] - vert[k]);
      }

      double prob_move_next = (1-beta);

      if(type[i] == vertex_type::max){
        model.addConstr(vert[i],GRB_GREATER_EQUAL,vert[j]*prob_move_next);
        model.addConstr(vert[i],GRB_GREATER_EQUAL,vert[k]*prob_move_next);
      }
      else if(type[i] == vertex_type::min){
        model.addConstr(vert[i],GRB_LESS_EQUAL,vert[j]*prob_move_next );
        model.addConstr(vert[i],GRB_LESS_EQUAL,vert[k]*prob_move_next );
      }
      else if(type[i] == vertex_type::ave){
        model.addConstr(vert[i],GRB_EQUAL,.5 *prob_move_next * vert[k] + .5 * prob_move_next  * vert[j]);
      }
      else if(type[i] == vertex_type::sink_max){
        model.addConstr(vert[i],GRB_EQUAL,1);
      }
      else if(type[i] == vertex_type::sink_min){
        model.addConstr(vert[i],GRB_EQUAL,0);
      }
    }

    // find coefficient for each vertex in the objective function.
    std::vector<double> coeffs(n, 1.0);

    // set objective function
    model.setObjective(objective_fn, GRB_MINIMIZE);

    // solve model
    model.optimize();

    // build probability vector
    for(int i = 0; i<n; i++){
      prob.push_back(vert[i].get(GRB_DoubleAttr_X));
    }
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch(...) {
    std::cout << "Exception during optimization" << std::endl;
  }
  return prob;
}

std::vector<bool> SSG::quadratic_program(){
    std::vector<bool> s(n,0);
    reconstruct_strategy(s, QP());
    return s;
}

std::vector<bool> SSG::quadratic_program(std::vector<bool> s){
    reconstruct_strategy(s, QP());
    return s;
}

std::vector<bool> SSG::hoffman_karp_min_LP(){
    std::vector<bool> s(n,0);
    return hoffman_karp_min_LP(s);
}
std::vector<bool> SSG::hoffman_karp_min_LP(std::vector<bool> s){
    bool strategys_updated;
    do{
        strategys_updated = false;
        auto probs = optimize_min_LP(s);
        strategys_updated |= reconstruct_strategy_min(s, probs);
        strategys_updated |= switch_max(s, probs);
    }while(strategys_updated);
    return s;
}

std::vector<bool> SSG::hoffman_karp_max_LP(){
    std::vector<bool> s(n,0);
    return hoffman_karp_max_LP(s);
}
std::vector<bool> SSG::hoffman_karp_max_LP(std::vector<bool> s){
    bool strategys_updated;
    do{
        strategys_updated = false;
        auto probs = optimize_max_LP(s);
        strategys_updated |= reconstruct_strategy_max(s, probs);
        strategys_updated |= switch_min(s, probs);
    }while(strategys_updated);
    return s;
}

std::vector<bool> SSG::hoffman_karp_LP(){
    std::vector<bool> s(n,0);
    return hoffman_karp_LP(s);
}
std::vector<bool> SSG::hoffman_karp_LP(std::vector<bool> s){
    if(min_vertices.size() < max_vertices.size()){
        return hoffman_karp_max_LP(s);
    }
    return hoffman_karp_min_LP(s);
}

std::vector<bool> SSG::hoffman_karp2(std::vector<bool> s){
    int iterations_to_convergence = 0;

    bool vert_switched_any;
    do{
        vert_switched_any = false;

        vert_switched_any |= optimize_min(s);
        auto probs = probabilities(s);
        
        iterations_to_convergence++;
        vert_switched_any = false;
        for(int cur_v: max_vertices){
            
            int cur_edge = outgoing_edge[cur_v][s[cur_v]];
            int other_edge = outgoing_edge[cur_v][!s[cur_v]];

            double p_cur = probs[cur_edge];
            double p_other = probs[other_edge];

            double delta = p_other - p_cur;
            delta = abs(delta);

            if(delta > tolerance && p_other > p_cur){
                // NOTICE: random()%2 has 1/2 of selecting current edge.
                vert_switched_any = true;
                s[cur_v] = !s[cur_v]; //switch edge.
            }
        }

    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}


std::vector<bool> SSG::hoffman_karp2_dermans(){
    std::vector<bool> s(n,0);
    return hoffman_karp2(s);
}
std::vector<bool> SSG::hoffman_karp2_dermans(std::vector<bool> s){
    int iterations_to_convergence = 0;

    optimize_min(s);
    optimize_max(s);

    bool vert_switched_any;
    do{
        vert_switched_any = false;

        vert_switched_any |= optimize_min(s);
        auto probs = probabilities(s);
        
        iterations_to_convergence++;
        vert_switched_any = false;
        for(int cur_v: max_vertices){
            
            int cur_edge = outgoing_edge[cur_v][s[cur_v]];
            int other_edge = outgoing_edge[cur_v][!s[cur_v]];

            double p_cur = probs[cur_edge];
            double p_other = probs[other_edge];

            double delta = p_other - p_cur;
            delta = abs(delta);

            if(delta > tolerance && p_other > p_cur){
                // NOTICE: random()%2 has 1/2 of selecting current edge.
                vert_switched_any = true;
                s[cur_v] = !s[cur_v]; //switch edge.
            }
        }

    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}

int SSG::optimize_min_iters(std::vector<bool> &s, std::vector<int> &count){
    auto p = probabilities(s);
    return optimize_min_iters(s,p,count);
}

int SSG::optimize_min_iters(std::vector<bool> &s, std::vector<double> probs, std::vector<int> &count){  
        bool min_switched_edge;
        //bool ever_switched = false;
        int min_verts_switched = 0;
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
                    count[cur_v] += 1;
                    probs = alternate_prob; //update probability vector
                    
                    min_switched_edge = true;
                    //ever_switched = true;
                    //std::cout << "n: " << s << std::endl;
                }
            }
            min_verts_switched+= min_vertices.size();
        }while(min_switched_edge);
    return min_verts_switched;
}

std::vector<bool> inverse_strategy(std::vector<bool> &s){
    for(int i = 0; i<s.size(); i++){
        s[i] = !s[i];
    }
    return s;
}


int SSG::hoffman_karp_n_iterations_inverse(){
    std::vector<bool> s = hoffman_karp();
    inverse_strategy(s);
    return hoffman_karp_n_iterations(s);
}


int SSG::hoffman_karp_n_iterations(){
    std::vector<bool> s(n,0);
    return hoffman_karp_n_iterations(s);
}
int SSG::hoffman_karp_n_iterations(std::vector<bool> s){
    int total_switches = 0;

    bool vert_switched_any;
    do{
        vert_switched_any = false;

        auto probs = probabilities(s);
        bool any_min_switched = optimize_min(s); //improve MIN
        vert_switched_any |= any_min_switched;

        probs = probabilities(s);
        //Improve MAX
        for(int cur_v: max_vertices){
            
            int cur_edge = outgoing_edge[cur_v][s[cur_v]];
            int other_edge = outgoing_edge[cur_v][!s[cur_v]];

            double p_cur = probs[cur_edge];
            double p_other = probs[other_edge];

            double delta = p_other - p_cur;
            delta = abs(delta);

            if(delta > tolerance && p_other > p_cur){
                    vert_switched_any = true;
                    s[cur_v] = !s[cur_v]; //switch edge.
            }
        }
        total_switches++;

    }while(vert_switched_any);

    return total_switches;
}

int SSG::hoffman_karp_n_iterations_dermans(){
    std::vector<bool> s(n,0);
    return hoffman_karp_n_iterations_dermans(s);
}
int SSG::hoffman_karp_n_iterations_dermans(std::vector<bool> s){

    optimize_min(s);
    optimize_max(s);

    int total_switches = 0;

    bool vert_switched_any;
    do{
        vert_switched_any = false;

        auto probs = probabilities(s);
        bool any_min_switched = optimize_min(s); //improve MIN
        vert_switched_any |= any_min_switched;

        probs = probabilities(s);
        //Improve MAX
        for(int cur_v: max_vertices){
            
            int cur_edge = outgoing_edge[cur_v][s[cur_v]];
            int other_edge = outgoing_edge[cur_v][!s[cur_v]];

            double p_cur = probs[cur_edge];
            double p_other = probs[other_edge];

            double delta = p_other - p_cur;
            delta = abs(delta);

            if(delta > tolerance && p_other > p_cur){
                    vert_switched_any = true;
                    s[cur_v] = !s[cur_v]; //switch edge.
            }
        }
        total_switches++;

    }while(vert_switched_any);

    return total_switches;
}


Eigen::MatrixXi compute_switch_matrix(std::vector<std::vector<double>> probs){
    if(!probs.size()){
        return Eigen::MatrixXi();
    }
    int n = probs[0].size();
    Eigen::MatrixXi times_changed = Eigen::MatrixXi::Zero(n-2, n-2);

    for(uint i = 1; i<probs.size(); i++){
         for(int r = 0; r<n-2; r++){
            for(int c = r+1; c<n-2; c++){
                double old_a = probs[i-1][c]; //index(a) > index(b)
                double old_b = probs[i-1][r];
                bool old_order = abs(old_a-old_b) > SSG::default_tolerance && old_a > old_b;

                double a = probs[i][c]; //index(a) > index(b)
                double b = probs[i][r];
                bool new_order = abs(a-b) > SSG::default_tolerance && a > b;
                if(new_order != old_order){
                    times_changed(r,c) += 1;
                    times_changed(c,r) += 1;
                }
            }
        }
    }
    return times_changed;
}



int SSG::hoffman_karp_n_iterations_print(){
    std::vector<bool> s(n,0);
    return hoffman_karp_n_iterations_print(s);
}
int SSG::hoffman_karp_n_iterations_print(std::vector<bool> s){
    int total_switches = 0;
    std::vector<int> switch_count(n,0);

    std::vector<std::vector<double>> n_probs(0); //before min makes improvements
    std::vector<std::vector<double>> x_probs(0); //before max makes improvesments
    
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        auto probs = probabilities(s);
        n_probs.push_back(probs);

        std::cout << "n: " << probs << " \tn: " << s << std::endl;
        int min_verts_switched = optimize_min_iters(s,switch_count); //improve MIN
        bool any_min_switched = min_verts_switched > (int)min_vertices.size();
        vert_switched_any |= any_min_switched;


        probs = probabilities(s);
        x_probs.push_back(probs);
        std::cout << "x: " << probs << " \tx: " << s << std::endl;
        //Improve MAX
        for(int cur_v: max_vertices){
            
            int cur_edge = outgoing_edge[cur_v][s[cur_v]];
            int other_edge = outgoing_edge[cur_v][!s[cur_v]];

            double p_cur = probs[cur_edge];
            double p_other = probs[other_edge];

            double delta = p_other - p_cur;
            delta = abs(delta);

            if(delta > tolerance && p_other > p_cur){
                    vert_switched_any = true;
                    s[cur_v] = !s[cur_v]; //switch edge.
                    switch_count[cur_v] += 1;
            }
        }
        total_switches++;

    }while(vert_switched_any);
    
    const char* sp = "      ";
    const char* type_c = " SNAXS";
    std:: cout << "v:  ";   for(int i = 0; i<n; i++){std::cout << i%10 << sp;} 
    std:: cout << " \tv: "; for(int i = 0; i<n; i++){std::cout << i%10 << " ";}                 std::cout << std::endl;
    std::cout << "s:  ";    for(int i = 0; i<n; i++){std::cout << switch_count[i] << sp;}
    std::cout << " \ts: ";  for(int i = 0; i<n; i++){std::cout << switch_count[i] << " ";}      std::cout << std::endl;
    std::cout << "t:  ";    for(int i = 0; i<n; i++){std::cout << type_c[(int)type[i]] << sp; } 
    std::cout << " \tt: ";  for(int i = 0; i<n; i++){std::cout << type_c[(int)type[i]] << " ";} std::cout << std::endl;
    std::cout << "\n\n";


    //keep track of how many times x&n switch each vertex.

    auto x_switches = compute_switch_matrix(x_probs);
    auto n_switches = compute_switch_matrix(n_probs);

    
    std::cout << "n:\n";
    std::cout << n_switches << "\n\n";
    std::cout << "x:\n";
    std::cout << x_switches << "\n\n";

    std::cout << switch_count<< std::endl;
    for(int i = 0; i<n; i++){std::cout << type_c[(int)type[i]] << " ";} std::cout << std::endl;
    for(int i = 0; i<n; i++){std::cout << i%10 << " ";}                 std::cout << std::endl;
    std::cout << std::endl;
    return total_switches;
}


std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> SSG::switch_matrices(){
    return switch_matrices(std::vector<bool>(n,0));
}
std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> SSG::switch_matrices(std::vector<bool> s){
    int total_switches = 0;
    std::vector<int> switch_count(n,0);

    std::vector<std::vector<double>> n_probs(0); //before min makes improvements
    std::vector<std::vector<double>> x_probs(0); //before max makes improvesments
    
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        auto probs = probabilities(s);
        n_probs.push_back(probs);

        std::cout << "n: " << probs << " \tn: " << s << std::endl;
        int min_verts_switched = optimize_min_iters(s,switch_count); //improve MIN
        bool any_min_switched = min_verts_switched > (int)min_vertices.size();
        vert_switched_any |= any_min_switched;


        probs = probabilities(s);
        x_probs.push_back(probs);
        std::cout << "x: " << probs << " \tx: " << s << std::endl;
        //Improve MAX
        for(int cur_v: max_vertices){
            
            int cur_edge = outgoing_edge[cur_v][s[cur_v]];
            int other_edge = outgoing_edge[cur_v][!s[cur_v]];

            double p_cur = probs[cur_edge];
            double p_other = probs[other_edge];

            double delta = p_other - p_cur;
            delta = abs(delta);

            if(delta > tolerance && p_other > p_cur){
                    vert_switched_any = true;
                    s[cur_v] = !s[cur_v]; //switch edge.
                    switch_count[cur_v] += 1;
            }
        }
        total_switches++;
    }while(vert_switched_any);

    auto x_switches = compute_switch_matrix(x_probs);
    auto n_switches = compute_switch_matrix(n_probs);

    return std::tuple<Eigen::MatrixXi, Eigen::MatrixXi>(n_switches, x_switches);
}


std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> SSG::final_switch_matrices(){
    return final_switch_matrices(std::vector<bool>(n,0));
}

std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> SSG::final_switch_matrices(std::vector<bool> s){
     int total_switches = 0;
    std::vector<int> switch_count(n,0);

    std::vector<std::vector<double>> n_probs(0); //before min makes improvements
    std::vector<std::vector<double>> x_probs(0); //before max makes improvesments
    
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        auto probs = probabilities(s);
        n_probs.push_back(probs);

        std::cout << "n: " << probs << " \tn: " << s << std::endl;
        int min_verts_switched = optimize_min_iters(s,switch_count); //improve MIN
        bool any_min_switched = min_verts_switched > (int)min_vertices.size();
        vert_switched_any |= any_min_switched;


        probs = probabilities(s);
        x_probs.push_back(probs);
        std::cout << "x: " << probs << " \tx: " << s << std::endl;
        //Improve MAX
        for(int cur_v: max_vertices){
            
            int cur_edge = outgoing_edge[cur_v][s[cur_v]];
            int other_edge = outgoing_edge[cur_v][!s[cur_v]];

            double p_cur = probs[cur_edge];
            double p_other = probs[other_edge];

            double delta = p_other - p_cur;
            delta = abs(delta);

            if(delta > tolerance && p_other > p_cur){
                    vert_switched_any = true;
                    s[cur_v] = !s[cur_v]; //switch edge.
                    switch_count[cur_v] += 1;
            }
        }
        total_switches++;
    }while(vert_switched_any);
    

    std::vector<std::vector<double>> last_x(x_probs.end()-2, x_probs.end());
    std::vector<std::vector<double>> last_n(n_probs.end()-2, n_probs.end());

    auto x_switches = compute_switch_matrix(last_x);
    auto n_switches = compute_switch_matrix(last_n);

    return std::tuple<Eigen::MatrixXi, Eigen::MatrixXi>(n_switches, x_switches);
}




std::vector<bool> SSG::incorrect_hoffman_karp(){
    std::vector<bool> s(n,0);
    return incorrect_hoffman_karp(s);
}

std::vector<bool> SSG::incorrect_hoffman_karp(std::vector<bool> s){
    auto probs = probabilities(s);

    bool vert_switched_any;

    std::vector<int> player_verts;
    do{
        vert_switched_any = false;

        vert_switched_any |= optimize_max(s,probs);
        vert_switched_any |= optimize_min(s);
        probs = probabilities(s);

    }while(vert_switched_any);
    return s;
}

std::vector<bool> SSG::tripathi_hoffman_karp(){
    std::vector<bool> s(n,0);
    return tripathi_hoffman_karp(s);  
}

std::vector<bool> SSG::tripathi_hoffman_karp(std::vector<bool> s){
    bool vert_switched_any;
    do{
        vert_switched_any = false;

        bool max_has_switchable_edge = false;
        bool max_vert_was_switched = false;
        auto probs = probabilities(s);
        do{ 
            vert_switched_any = false;
            for(int cur_v: max_vertices){
                
                int cur_edge = outgoing_edge[cur_v][s[cur_v]];
                int other_edge = outgoing_edge[cur_v][!s[cur_v]];

                double p_cur = probs[cur_edge];
                double p_other = probs[other_edge];

                double delta = p_other - p_cur;
                delta = abs(delta);

                if(delta > tolerance && p_other > p_cur){
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
        
        vert_switched_any |= optimize_min(s);

        probs = probabilities(s);
        
    }while(vert_switched_any);
    //std::cout << "its: " << iterations_to_convergence << "  ave player loops: " << player_loops / (double)(2*iterations_to_convergence) << std::endl;
    return s;
}

std::vector<bool> SSG::tripathi_hoffman_karp_LP(){
    std::vector<bool> s(n,0);
    return tripathi_hoffman_karp_LP(s);  
}

std::vector<bool> SSG::tripathi_hoffman_karp_LP(std::vector<bool> s){
    bool strategys_updated;
    do{
        strategys_updated = false;
        auto probs = optimize_min_LP(s);
        strategys_updated |= reconstruct_strategy_min(s, probs);
        strategys_updated |= switch_max_tripathi(s, probs);
    }while(strategys_updated);
    return s;
}

std::vector<bool> SSG::ludwig(){
    std::vector<bool> s(n,0);
    return ludwig(s);
}

std::vector<bool> SSG::ludwig(const std::vector<bool> &s){
    std::vector<bool> selected(n,0);
    return s;
}

std::vector<bool> SSG::ludwig_iterative_LP(){
    std::vector<bool> s(n,0);
    return ludwig_iterative_LP(s);
}

std::vector<bool> SSG::ludwig_iterative_LP(std::vector<bool> s){
    std::vector<int> s_removed(max_vertices.begin(), max_vertices.end());
    std::vector<int> s_in_game(0);

    std::shuffle(std::begin(s_removed), std::end(s_removed), rng);

    auto p = optimize_min_LP(s);
    reconstruct_strategy(s, p);

    while(s_removed.size()){

        int cur_s = s_removed.back(); s_removed.pop_back();

        s_in_game.push_back(cur_s);

        for(uint k = 0; k<s_in_game.size(); k++){
            int j = s_in_game[k];

            s[j] = !s[j];
            auto new_p = probabilities(s);
            s[j] = !s[j];

            double delta = abs(p[j] - new_p[j]);

            if(delta > tolerance && new_p[j] > p[j]){
                s[cur_s] = !s[cur_s];
                p = optimize_min_LP(s);
                reconstruct_strategy(s, p);

                std::shuffle(std::begin(s_in_game), std::end(s_in_game),rng);

                while(s_in_game.size()>0){
                    int v = s_in_game.back();
                    s_in_game.pop_back();
                    s_removed.push_back(v);
                }
                break;
            }
        }
    } //std::cout << std::endl;

    return s;
}




std::vector<bool> SSG::ludwig_iterative(){
    std::vector<bool> s(n,0);
    return ludwig_iterative(s);
}

std::vector<bool> SSG::ludwig_iterative(std::vector<bool> s){
    std::vector<int> s_removed(max_vertices.begin(), max_vertices.end());
    std::vector<int> s_in_game(0);

    std::shuffle(std::begin(s_removed), std::end(s_removed), rng);

    optimize_min(s);

    auto p = probabilities(s);

    while(s_removed.size()){

        int cur_s = s_removed.back(); s_removed.pop_back();
        
        s_in_game.push_back(cur_s);


            for(uint k = 0; k<s_in_game.size(); k++){
            int j = s_in_game[k];

            s[j] = !s[j];
            auto new_p = probabilities(s);
            s[j] = !s[j];

            double delta = abs(p[j] - new_p[j]);

            if(delta > tolerance && new_p[j] > p[j]){
                s[cur_s] = !s[cur_s];
                optimize_min(s);
                p = probabilities(s);
                
                std::shuffle(std::begin(s_in_game), std::end(s_in_game),rng);

                while(s_in_game.size()>0){
                    int v = s_in_game.back();
                    s_in_game.pop_back();
                    s_removed.push_back(v);
                }
                break;
            }
            }
    } //std::cout << std::endl;

    return s;
}

std::vector<bool> SSG::ludwig_iterative_OG(){
    std::vector<bool> s(n,0);
    return ludwig_iterative_OG(s);
}

std::vector<bool> SSG::ludwig_iterative_OG(std::vector<bool> s){
    std::vector<int> s_removed(max_vertices.begin(), max_vertices.end());
    std::vector<int> s_in_game(0);

    std::shuffle(std::begin(s_removed), std::end(s_removed), rng);

    optimize_min(s);

    while(s_removed.size()){

        int cur_s = s_removed.back(); s_removed.pop_back();

        s_in_game.push_back(cur_s);

        auto p = probabilities(s);

        for(uint k = 0; k<s_in_game.size(); k++){
            int j = s_in_game[k];

            s[j] = !s[j];
            auto new_p = probabilities(s);
            s[j] = !s[j];

            double delta = abs(p[j] - new_p[j]);

            if(delta > tolerance && new_p[j] > p[j]){
                s[cur_s] = !s[cur_s];
                optimize_min(s);

                std::shuffle(std::begin(s_in_game), std::end(s_in_game),rng);

                while(s_in_game.size()>0){
                    int v = s_in_game.back();
                    s_in_game.pop_back();
                    s_removed.push_back(v);
                }
                break;
            }
        } //std::cout << std::endl;
    } //std::cout << std::endl;

    return s;
}




SSG SSG::random_game_loopless(int n){
    SSG game(n);

    if(n<2)
        return game;

    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);

    const std::vector<vertex_type> v_type= {vertex_type::min, vertex_type::ave, vertex_type::max};

    for(int i = n-3; i>= 0; i--){
        vertex_type type = v_type[rand()%3];
        int j = random()%(n - i) + i;
        int k = random()%(n - i) + i;

        game.set_vertex(i, type, j, k);
    }

    return game;
}

SSG SSG::random_game(int n){
    SSG game(n);

    if(n<2)
        return game;

    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);

    const std::vector<vertex_type> v_type= {vertex_type::min, vertex_type::ave, vertex_type::max};

    for(int i = n-3; i>= 0; i--){
        vertex_type type = v_type[rand()%3];
        int j = random()%n;
        int k = random()%n;

        game.set_vertex(i, type, j, k);
    }

    return game;
}

SSG SSG::random_nontrivial_game(int n){
    SSG game(n);

    if(n<2)
        return game;

    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);

    const std::vector<vertex_type> v_type= {vertex_type::min, vertex_type::ave, vertex_type::max};

    std::vector<int> nax(3,0); //find how many min, ave, max verts there will be.
    for(int i = 0; i<n-2; i++){
        nax[random()%3]++;
    }

    nax[1] += nax[0];
    nax[2] += nax[1];

    int k = 0;

    for(int v = 0; v< n-2; v++){

        while(v >= nax[k]){
            k++;
        }

        vertex_type type = v_type[k];

        int e1;
        int e2;

        std::set<int> exclude;
        exclude.insert(v);

        switch(type){
            case ave:
            {
                for(int i_v: game.incoming_edges[v]){
                    if(game.type[i_v] == vertex_type::ave){
                        exclude.insert(i_v);
                    }
                }
                e1 = random_in_range_exclude(0, n, exclude);
                exclude.insert(e1);
                e2 = random_in_range_exclude(0, n, exclude);
                break;
            }
            case min:
            {
                for(int i_v: game.incoming_edges[v]){
                    if(game.type[i_v] == vertex_type::min || game.type[i_v] == vertex_type::max){
                        exclude.insert(i_v);
                    }
                }
                e1 = random_in_range_exclude(0, n-2, exclude);
                exclude.insert(e1);
                e2 = random_in_range_exclude(0, n-2, exclude);
                break;
            }
            case max:
            {
                for(int i_v: game.incoming_edges[v]){
                    if(game.type[i_v] == vertex_type::max || game.type[i_v] == vertex_type::min){
                        exclude.insert(i_v);
                    }
                }
                e1 = random_in_range_exclude(0, n-2, exclude);
                exclude.insert(e1);
                e2 = random_in_range_exclude(0, n-2, exclude);
                break;
            }
        }
        game.set_vertex(v, type, e1, e2);
    }

    return game;
}

SSG SSG::random_game_equal_split(int n){
    SSG game(n);

    if(n<2)
        return game;

    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);

    int p_verts = (n-2)/3;

    for(int i = n-3; i> n-3 - p_verts; i--){
        vertex_type type = vertex_type::max;
        int j = random()%n;
        int k = random()%n;

        game.set_vertex(i, type, j, k);
    }

    for(int i = n-3 - p_verts; i> n-3 - (p_verts * 2); i--){
        vertex_type type = vertex_type::ave;
        int j = random()%n;
        int k = random()%n;

        game.set_vertex(i, type, j, k);
    }
    for(int i = n-3 - (p_verts*2); i>= 0; i--){
        vertex_type type = vertex_type::min;
        int j = random()%n;
        int k = random()%n;

        game.set_vertex(i, type, j, k);
    }

    return game;
}

SSG SSG::random_game_n_max(int n_max, int n){
    SSG game(n);

    if(n<2)
        return game;

    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1); 

    const std::vector<vertex_type> v_type= {vertex_type::min, vertex_type::ave, vertex_type::max};
    for(int i = n-3; i>=0; i--){
        vertex_type type = v_type[random()%2];
        int j = random()%n;
        int k = random()%n;

        game.set_vertex(i, type, j, k);
    }

    auto maxes = permutation_in_range(0,n-2);

    for(int i = 0; i<n_max; i++){
        int cur_max = maxes.back(); maxes.pop_back();

        int j = random()%n;
        int k = random()%n;

        game.force_vertex(cur_max, vertex_type::max, j, k);
    }

    return game;
}


SSG SSG::double_chain(int n){
    SSG game(n);

    if(n<2)
        return game;

    int last_min_vertex = (n-2)/2;


    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);

    game.set_vertex(n-3, vertex_type::max, n-2, n-1);

    for(int i = n-4; i>=last_min_vertex+1; i--){
        vertex_type type = vertex_type::max;
        int j = n-2;
        int k = i+1;

        game.set_vertex(i, type, j, k);
    }

    game.set_vertex(last_min_vertex, vertex_type::min, n-1, n-2);

    for(int i = last_min_vertex-1; i>=0; i--){
        vertex_type type = vertex_type::min;
        int j = n-1;
        int k = i+1;

        game.set_vertex(i, type, j, k);
    }

    return game;
}

SSG SSG::hard_game_max(int n){
    SSG game(n);

    if(n<2)
        return game;

    game.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);
    game.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);

    game.set_vertex(n-3, vertex_type::max, 0, n-1);

    for(int i = n-4; i>=0; i--){
        vertex_type type = vertex_type::max;
        int j = 0;
        int k = i+1;

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


SSG SSG::read_game(std::string file_name){
    std::ifstream file;
    file.open(file_name);
    SSG g = SSG::read_game_file(file);
    file.close();
    return g;
}





SSG SSG::read_game_expand(std::string file_name, int extra_vertices){
    std::ifstream file;
    file.open(file_name);
    SSG g = SSG::read_game_file_expand(file, extra_vertices);
    file.close();
    return g;
}


SSG SSG::read_game_file_expand(std::ifstream &file, int extra_vertices){
    std::string cur_line;

    std::getline(file, cur_line);

    int n = stoi(cur_line);
    SSG game(n+extra_vertices);

    for(int i = 0; i<n; i++){
        std::getline(file, cur_line);

        int v, type, e1, e2;
        std::istringstream iss(cur_line);
        iss >> v >> type >> e1 >> e2;

        #define v_macro(x) x = x>n-3? x+extra_vertices : x

        v_macro(v);
        v_macro(e1);
        v_macro(e2);

        game.set_vertex(v,(vertex_type)type,e1,e2);
    }
    return game;
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


bool SSG::probs_match(const std::vector<double> &p1, const std::vector<double> &p2, double tolerance){
    for(uint i = 0; i< std::min(p1.size(), p2.size());i++){
        double delta = std::abs(p1[i] - p2[i]);
        if(delta > tolerance){
            return false;
        }
    }
    return true;
}

// increments the game vertices. used for bruteforce game generation
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

    type = type_overflow ? vertex_type::min : type;

    //cout << v << " -> "<< e1 << " " << e2 << endl;
    g.force_vertex(v, type, e1, e2);

    
    if(!type_overflow && type == vertex_type::max){
        if(e1 == e2)
            return next_combination_most(g,v);
    }

    if(!type_overflow && type == vertex_type::min){
        if(e1 == e2 || e1 == v || e2 == v)
            return next_combination_most(g,v);
    }
    

    return type_overflow ? next_combination_most(g, v-1) : true;
}

//bruteforces hardest game
SSG SSG::hardest_possible_game(int n){
    SSG g(n);
    for(int i = 0; i<n-2; i++){
        g.set_vertex(i,vertex_type::ave,0,0);
    }
    g.set_vertex(n-1, vertex_type::sink_max, n-1, n-1);
    g.set_vertex(n-2, vertex_type::sink_min, n-2, n-2);

    double max = 0;
    std::vector<SSG> max_gs(0, SSG(n));

    do{
        if(g.max_vertices.size()>0){
            double cur = g.hoffman_karp_n_iterations();

            if(cur/(double)g.max_vertices.size() > max){
                max = cur;
                max_gs.push_back(g.copy());
            }
        }
    }while(next_combination_most(g, n-3));

    return max_gs.back();
}



std::ostream& operator<<(std::ostream& os, const SSG &game)
{
    os << game.n << std::endl;
    for(int i = 0; i<game.n; i++){
        os << i << "\t" << game.type[i] << "\t" << game.outgoing_edge[i][0] << "\t" <<game.outgoing_edge[i][1] << std::endl;
    }
    return os;
}

// << operator for strategies
std::ostream& operator<<=(std::ostream& os, const SSG &game){
    os << game.n << std::endl;
    for(int i = 0; i<game.n-2; i++){
        os << i <<"_"<< vertex_type_names[game.type[i]] << "\t" << game.outgoing_edge[i][0] << "_" << vertex_type_names[game.type[game.outgoing_edge[i][0]]] << std::endl;
        os << i <<"_"<< vertex_type_names[game.type[i]] << "\t" << game.outgoing_edge[i][1] << "_" << vertex_type_names[game.type[game.outgoing_edge[i][1]]] << std::endl;
    }

    os << game.n-2 << "_" << vertex_type_names[game.type[game.n-2]] << std::endl;
    os << game.n-1 << "_" << vertex_type_names[game.type[game.n-1]];

    return os;
}


//PRIVATE SETTERS
void SSG::set_edges(int vertex, int e1, int e2){

    /*if(e1 > e2){
        int e3 = e1;
        e1 = e2;
        e2 = e3;
    }*/

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

    incoming_edges[e1].insert(vertex);
    incoming_edges[e2].insert(vertex);
}

void SSG::force_edges(int vertex, int e1, int e2){
    if(e1 > e2){
        int e3 = e1;
        e1 = e2;
        e2 = e3;
    }

    //out of range for vertices
    if(e2 >= n || e1 < 0 || e1 >= n || e2 < 0){
        std::cerr << "WARNING SSG::set_edges e1 or e2 out of bounds. e1: "  << e1 << "\te2: " << e2 << std::endl;
        return;
    }

    int old_e1 = outgoing_edge[vertex][0];
    int old_e2 = outgoing_edge[vertex][0];

    incoming_edges[old_e1].erase(vertex);
    incoming_edges[old_e2].erase(vertex);

    outgoing_edge[vertex][0] = e1;
    outgoing_edge[vertex][1] = e2;

    incoming_edges[e1].insert(vertex);
    incoming_edges[e2].insert(vertex);
}

void SSG::set_vertex(int vertex, vertex_type type, int e1, int e2){
    (*this).set_vertex_type(vertex, type);
    (*this).set_edges(vertex, e1, e2);
}

void SSG::force_vertex(int vertex, vertex_type type, int e1, int e2){
    (*this).force_vertex_type(vertex, type);
    (*this).force_edges(vertex, e1, e2);
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

std::ostream& operator<<(std::ostream& os, const std::vector<int> &vec){
    for(auto it = vec.begin(); it!=vec.end(); it++){
        os << *it << " ";
    }
    return os;
}


void SSG::print_vertex_types(){
    for(int i = 0; i<n; i++){
        switch(type[i]){
            std::cout << i << ":";
            case vertex_type::min:{
                std::cout << "n";
                break;
            }
            case vertex_type::ave:{
                std::cout << "a"; 
                break;
            }
            case vertex_type::max:{
                std::cout << "x"; 
                break;
            }
            case vertex_type::sink_min:{
                std::cout << "N"; 
                break;
            }
            case vertex_type::sink_max:{
                std::cout << "X"; 
                break;
            }
        }
        std::cout << " ";
    }
    std::cout << std::endl;
}