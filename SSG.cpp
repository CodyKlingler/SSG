#include <iostream>
#include <random>
#include "SSG.h"
#include "Strategy.h"

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
    int min_sink_vertex = -1;
    int max_sink_vertex = -1;

    strat = NULL;
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
    //out of range for vertices
    if(e1 >= n || e1 < 0 && e2 >= n || e1 < 0){
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


void SSG::start(int starting_vertex, Strategy combined_strategy){
    token = starting_vertex;
    strat = combined_strategy;
    n_steps = 0;
}

int SSG::step(int n_steps = 1){
    for(int i = 0; i<n_steps; i++){
        switch(type[token]){
            case undefined:
                std::cerr << "WARNING: SSG::step token at vertex " << token << " is undefined" << std::endl;
                return -999;
            case ave:
                int next_edge = random()%2;
                token = outgoing_edge[token][next_edge];
                n_steps++;
                break;
            case sink_max:
                return 1;
            case sink_min:
                return 0;
            default:
              token = strat[token];
              n_steps++;
              break;
        }

        if(n_steps > 1000*n){
            return -1;
        }
    }
    return -1;
}

int SSG::play(int starting_vertex, Strategy combined_strategy){
    start(starting_vertex, combined_strategy);
}

double SSG::play_n(int starting_vertex, Strategy combined_strategy, int n_trials){

}

void SSG::print_graph(){
    
    std::vector<std::vector<int>> vec_list = {min_vertices, max_vertices, ave_vertices};
    std::vector<const char*> vec_names =     {"min_vertices", "max_vertices", "ave_vertices"};

    for(int i = 0; i< vec_list.size(); i++){

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
