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

//creates SSG without loops
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

std::vector<bool> SSG::hoffman_karp(){
    std::vector<bool> s(n,0);
    return hoffman_karp(s);
}

//checks if u's prob depends on v's
bool SSG::probability_depends_on(int u, int v){
    return false;
}


int xxx = 0;
const int nnn = 10000000;

bool print_prob = false;

bool exported_b4 = false;

std::vector<bool> SSG::hoffman_karp(std::vector<bool> s){
    
    auto probs = probabilities(s);

    bool vert_switched_any;

    std::vector<int> player_verts;
    
    do{ // while(vert_switched_any);
        vert_switched_any = false;
        if(xxx > nnn) 
            std::cout << "&&" << std::endl;
        player_verts = (player_verts == min_vertices)?max_vertices:min_vertices; //switch player vector
        std::vector<std::vector<int>> player_verts_vec = {min_vertices, max_vertices};
        for(auto player_verts: player_verts_vec){ 
            if(xxx > nnn) 
                std::cout << "||" << std::endl;
            for(auto it = player_verts.begin(); it!= player_verts.end(); it++){
                int cur_v = *it;
                int cur_edge = outgoing_edge[cur_v][s[cur_v]];
                int other_edge = outgoing_edge[cur_v][!s[cur_v]];
                
                s[cur_v] = !s[cur_v];
                auto alternate_prob = probabilities(s);
                s[cur_v] = !s[cur_v];

                double p_cur = probs[cur_v];
                double p_other = alternate_prob[cur_v];

                /*
                if(p_cur > 1.05){
                    std::cout << "probs: " << probs << std::endl;
                    print_prob = true;
                    probabilities(s);
                    print_prob = false;
                }

                if(p_other > 1.05){
                    std::cout << "other probs: " << alternate_prob << std::endl;
                    print_prob = true;
                    s[cur_v] = !s[cur_v];
                    auto alternate_prob = probabilities(s);
                    s[cur_v] = !s[cur_v];
                    print_prob = false;
                }
                */


                double delta = p_other - p_cur;
                delta = abs(delta);
                if(xxx++ > nnn){

                    if(!exported_b4){
                        exported_b4 = true;

                        std::ofstream myfile;
                        myfile.open("this.txt");

                        myfile << *this;
                        myfile.close();
                    }

                    printf("cur_v: %i\t ",cur_v);
                    printf(vertex_type_names[type[cur_v]]);
                    printf("\t cur_p(%i): %f\t other_p(%i): %f delta: %f\n", cur_edge, p_cur, other_edge, p_other, delta);
                    std::cout << "cur prob arr: " << probs << std::endl;
                    std::cout << "alt prob arr: " << alternate_prob << std::endl;
                    std::cout << "strategy" << s << std::endl;
                }
                    
                if(delta > .001){

                    if(p_other > p_cur && player_verts == max_vertices){
                        vert_switched_any = true;
                        s[cur_v] = !s[cur_v]; //switch edge.
                        probs = alternate_prob; //update probability vector
                        it = player_verts.begin();
                        if(xxx++ > nnn){
                            std::cout << "swap max" << std::endl;
                        }
                    }
                    else if (p_other < p_cur && player_verts == min_vertices){
                        vert_switched_any = true;
                        s[cur_v] = !s[cur_v]; //switch edge.
                        probs = alternate_prob; //update probability vector
                        it = player_verts.begin();
                        if(xxx++ > nnn)
                            std::cout << "swap min" << std::endl;
                    }
                }
            }
        }
    }while(vert_switched_any);
    if(xxx++ > nnn)
        std::cout << "EXIT" << std::endl;
    return s;
}




std::vector<double> SSG::probabilities(std::vector<bool> combined_strategy){
    using namespace Eigen;

    std::vector<Triplet<double>> coeffs(n*3);

    for(std::vector<int>::iterator it = max_vertices.begin(); it != max_vertices.end(); it++){
        //mat.coeffRef(*it,*it) = 1;
        bool b = combined_strategy[*it];
        int p = outgoing_edge[*it][b];
        //mat.coeffRef(*it,p) = -1;

        coeffs.push_back(Triplet<double>(*it,*it,1));
        if(*it != p){
            coeffs.push_back(Triplet<double>(*it,p,-1));
        }
    }

    for(std::vector<int>::iterator it = min_vertices.begin(); it != min_vertices.end(); it++){
        //mat.coeffRef(*it,*it) = 1;
        bool b = combined_strategy[*it];
        int p = outgoing_edge[*it][b];
        //mat.coeffRef(*it,p) = -1;

        coeffs.push_back(Triplet<double>(*it,*it,1));
        if(*it != p){
            coeffs.push_back(Triplet<double>(*it,p,-1));
        }
    }

    for(std::vector<int>::iterator it = ave_vertices.begin(); it != ave_vertices.end(); it++){
        //mat.coeffRef(*it,*it) = 1;

        int p1 = outgoing_edge[*it][0];
        int p2 = outgoing_edge[*it][1];
        //mat.coeffRef(*it,p1) = -.5;
        //mat.coeffRef(*it,p2) = -.5;

        coeffs.push_back(Triplet<double>(*it,*it,1));
        if(*it != p1){
            coeffs.push_back(Triplet<double>(*it,p1,-.5));
        }
        if(*it != p2 && p2 != p1){
            coeffs.push_back(Triplet<double>(*it,p2,-.5));
        }
        
    }

    

    //mat.coeffRef(min_sink_vertex,min_sink_vertex) = 1;
    //mat.coeffRef(max_sink_vertex, max_sink_vertex) = 1;



    coeffs.push_back(Triplet<double>(max_sink_vertex,max_sink_vertex,1));

    coeffs.push_back(Triplet<double>(min_sink_vertex,min_sink_vertex,1));
    
    SparseMatrix<double> mat(n,n); 
    //mat.reserve(MatrixXi::Constant(n,3));
    mat.setFromTriplets(coeffs.begin(), coeffs.end());  // fill A and b;
    //BiCGSTAB<SparseMatrix<double> > solver;

    LeastSquaresConjugateGradient<SparseMatrix<double> > solver;
    mat.makeCompressed();

    
    
    solver.compute(mat);

    VectorXd vec = VectorXd::Zero(n);
    vec(max_sink_vertex) = 1;

    //Use the factors to solve the linear system 
    VectorXd x = solver.solve(vec).cwiseAbs(); 

    if(0 && print_prob){
        std::cout << mat << std::endl;
        std::cout << "linear system solver:   estimated error: " << solver.error() << std::endl;
    }
    
    double* probs_temp = x.data();

    int vec_size = x.rows();
 
    return std::vector<double>(probs_temp, probs_temp + vec_size);
}



SSG::SSG(int vertices){
    n = vertices;
    type = new vertex_type[vertices]{undefined};

    outgoing_edge = (int**) malloc(sizeof(int**)*vertices);
    for(int i = 0; i<vertices; i++){
        outgoing_edge[i] = (int*)malloc(sizeof(int) * 2);
        outgoing_edge[i][0] = -1;
        outgoing_edge[i][1] = -1;
    }

    incoming_edges = new std::vector<int>[vertices];

    min_vertices = *new std::vector<int>;
    max_vertices = *new std::vector<int>;
    ave_vertices = *new std::vector<int>;
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


void SSG::start(int starting_vertex, std::vector<bool> combined_strategy){
    token = starting_vertex;
    strat = combined_strategy;
    n_steps = 0;
}

/*
    1: max wins
    0: min wins
    -1: no winner yet

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


std::ostream& operator<<(std::ostream& os, const SSG &game)
{
    os << game.n << std::endl;
    for(int i = 0; i<game.n; i++){
        os << i << "\t" << game.type[i] << "\t" << game.outgoing_edge[i][0] << "\t" <<game.outgoing_edge[i][1] << std::endl;
    }
    return os;
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

// << operator for strategies
std::ostream& operator<<(std::ostream& os, const std::vector<bool> &vec){
    for(auto it = vec.begin(); it!=vec.end(); it++){
        os << *it << " ";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<double> &vec){
    for(auto it = vec.begin(); it!=vec.end(); it++){
        os << *it << " ";
    }
    return os;
}
