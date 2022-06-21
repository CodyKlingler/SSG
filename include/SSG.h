#pragma once 

#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iostream>

#include "Strategy.h"

enum vertex_type{
    undefined = 0,
    ave = 1,
    min = 2,
    max = 3,
    sink_min = 4,
    sink_max = 5
};

//const char** vertex_type_names;


class SSG{
    public:
        int n;

        std::vector<std::vector<int>> outgoing_edge;
        std::vector<std::vector<int>> incoming_edges;
        std::vector<vertex_type> type;
        std::vector<int> min_vertices;
        std::vector<int> max_vertices;
        std::vector<int> ave_vertices;
        int min_sink_vertex = -1;
        int max_sink_vertex = -1;
        int token;
        std::vector<bool> strat;
        double tolerance = .01;
        const int c = 2;
        double beta = 0.00010;
   
        SSG(int n_vertices);
        //SSG(const SSG &game); //TODO
        //SSG(std::ifstream& file); //TODO
        void set_vertex(int vertex, vertex_type type, int e1, int e2);


        void start(int starting_vertex, std::vector<bool> strategy);
        int step(int n_steps = 1);
        int play(int starting_vertex, std::vector<bool> strategy);
        double play_n(int starting_vertex, std::vector<bool> strategy, int n_trials);

        void print_graph();

        SSG stopping_game();

        std::vector<double> probabilities(const std::vector<bool> &strategy);
        std::vector<double> exact_probabilities(const std::vector<bool> &strategy);


        //STRATEGY SOLVING ALGORITHMS

        bool optimize_min(std::vector<bool> &s);
        bool optimize_min(std::vector<bool> &s, std::vector<double> probs);

        bool optimize_max(std::vector<bool> &s);
        bool optimize_max(std::vector<bool> &s, std::vector<double> probs);

        std::vector<bool> bruteforce();
        std::vector<bool> bruteforce(std::vector<bool> s);

        std::vector<bool> hoffman_karp();
        std::vector<bool> hoffman_karp(std::vector<bool> s);
        
        std::vector<bool> ludwig();
        std::vector<bool> ludwig(const std::vector<bool> &s);

        std::vector<bool> ludwig_iterative();
        std::vector<bool> ludwig_iterative(std::vector<bool> s);

        std::vector<bool> incorrect_hoffman_karp(std::vector<bool> s);
        std::vector<bool> incorrect_hoffman_karp();

        std::vector<bool> tripathi_hoffman_karp(std::vector<bool> s);
        std::vector<bool> tripathi_hoffman_karp();

        //GAME AND STRATEGY GENERATION
        static SSG random_game_loopless(int n);
        static SSG random_game(int n);
        static std::vector<bool> random_strategy(int n);
        static SSG read_game_file(std::ifstream &file);
        static std::vector<bool> read_strategy_file(std::ifstream &file);
        static bool probs_match(const std::vector<double> &p1, const std::vector<double> &p2, double tolerance);


        friend std::ostream& operator<<(std::ostream& stream, const SSG& game);

    private:
        int n_steps;
        int n_steps_terminate; //after this many steps, the game will be won by MIN.

        void set_vertex_type(int vertex, vertex_type type);
        void set_edges(int vertex, int e1, int e2);
        bool increment_min_strat(std::vector<bool> &s);
        bool increment_max_strat(std::vector<bool> &s);
        std::vector<bool> bruteforce_max(std::vector<bool> s);
};

// << operator for strategies
std::ostream& operator<<(std::ostream& os, const std::vector<bool> &vec);

// << operator for strategies
std::ostream& operator<<(std::ostream& os, const std::vector<double> &vec);