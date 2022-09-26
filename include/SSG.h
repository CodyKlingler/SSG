#pragma once 

#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <tuple>
#include <set>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Strategy.h"

enum vertex_type{
    undefined = 0,
    sink_min = 1,
    min = 2,
    ave = 3,
    max = 4,
    sink_max = 5
};

//const char** vertex_type_names;


class SSG{
    public:
        int n;

        std::vector<std::vector<int>> outgoing_edge;
        std::vector<std::set<int>> incoming_edges;
        std::vector<vertex_type> type;
        std::vector<int> min_vertices;
        std::vector<int> max_vertices;
        std::vector<int> ave_vertices;
        int min_sink_vertex = -1;
        int max_sink_vertex = -1;
        int token;
        std::vector<bool> strat;
        constexpr static double default_tolerance = .0075;
        double tolerance = default_tolerance;
        const int c = 2;
        double beta = 0.0010;
        

        SSG(int n_vertices);
        //SSG(const SSG &game); //TODO
        //SSG(std::ifstream& file); //TODO
        void set_vertex(int vertex, vertex_type type, int e1, int e2);
        void force_vertex(int vertex, vertex_type type, int e1, int e2);

        void start(int starting_vertex, std::vector<bool> strategy);
        int step(int n_steps = 1);
        int play(int starting_vertex, std::vector<bool> strategy);
        double play_n(int starting_vertex, std::vector<bool> strategy, int n_trials);

        void print_graph();

        SSG copy();

        SSG reduce();

        SSG stopping_game();

        std::vector<double> probabilities(const std::vector<bool> &strategy);
        std::vector<double> exact_probabilities(const std::vector<bool> &strategy);


        //STRATEGY SOLVING ALGORITHMS

        //reconstruct min from prob vector. return true if any changes made to strategy
        bool reconstruct_strategy(std::vector<bool> &strategy, std::vector<double> p);
        bool switch_max(std::vector<bool> &strategy, std::vector<double> p);
        std::vector<double> optimize_min_LP(std::vector<bool> &strategy);


        bool optimize_min(std::vector<bool> &s);
        bool optimize_min(std::vector<bool> &s, std::vector<double> probs);

        int optimize_min_iters(std::vector<bool> &s,std::vector<int> &count);
        int optimize_min_iters(std::vector<bool> &s, std::vector<double> probs,std::vector<int> &count);

        bool optimize_max(std::vector<bool> &s);
        bool optimize_max(std::vector<bool> &s, std::vector<double> probs);

        std::vector<bool> strategy_prefer_own();

        std::vector<bool> bruteforce();
        std::vector<bool> bruteforce(std::vector<bool> s);

        std::vector<bool> hoffman_karp();
        std::vector<bool> hoffman_karp(std::vector<bool> s);
        std::vector<bool> hoffman_karp2();
        std::vector<bool> hoffman_karp2(std::vector<bool> s); //does min before max.
        std::vector<bool> hoffman_karp2_dermans();
        std::vector<bool> hoffman_karp2_dermans(std::vector<bool> s);
        std::vector<bool> hoffman_karp_LP();
        std::vector<bool> hoffman_karp_LP(std::vector<bool> s);


        int hoffman_karp_n_iterations();
        int hoffman_karp_n_iterations(std::vector<bool> s);
        int hoffman_karp_n_iterations_inverse();
        int hoffman_karp_n_iterations_dermans();
        int hoffman_karp_n_iterations_dermans(std::vector<bool> s); //finds opt for min and max before doing hoffman karp.
        int hoffman_karp_n_iterations_print();
        int hoffman_karp_n_iterations_print(std::vector<bool> s);
        
        std::vector<bool> ludwig();
        std::vector<bool> ludwig(const std::vector<bool> &s);

        std::vector<bool> ludwig_iterative();
        std::vector<bool> ludwig_iterative(std::vector<bool> s);

        std::vector<bool> ludwig_iterative_2();
        std::vector<bool> ludwig_iterative_2(std::vector<bool> s);

        std::vector<bool> incorrect_hoffman_karp(std::vector<bool> s);
        std::vector<bool> incorrect_hoffman_karp();

        std::vector<bool> tripathi_hoffman_karp(std::vector<bool> s);
        std::vector<bool> tripathi_hoffman_karp();


        //[0] is min matrix, [1] is max matrix
        std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> switch_matrices();
        std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> switch_matrices(std::vector<bool> s);
        std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> final_switch_matrices();
        std::tuple<Eigen::MatrixXi, Eigen::MatrixXi> final_switch_matrices(std::vector<bool> s);

        //GAME AND STRATEGY GENERATION
        static SSG random_game_loopless(int n);
        static SSG random_game(int n);
        static SSG random_nontrivial_game(int n); // random game excluding some poly-time solvable sub-structures (loops of min-vertices, player vertices that point to sinks, self-edges, etc.)
        static SSG random_game_equal_split(int n);
        static SSG random_game_n_max(int n_max, int n);
        static SSG hardest_possible_game(int n);    //bruteforce the game that takes the highest number of iterations by hk_iters
        static SSG hardest_game_nontrivial(int n);
        static SSG hard_game_max(int n);

        static std::vector<bool> random_strategy(int n);
        static SSG read_game_file(std::ifstream &file);
        static SSG read_game_file_expand(std::ifstream &file, int extra_vertices);
        static SSG read_game(std::string file_name);
        static SSG read_game_expand(std::string file_name, int extra_vertices);
        static std::vector<bool> read_strategy_file(std::ifstream &file);
        static bool probs_match(const std::vector<double> &p1, const std::vector<double> &p2, double tolerance);

        friend std::ostream& operator<<(std::ostream& stream, const SSG& game);
        friend std::ostream& operator<<=(std::ostream& stream, const SSG& game);

    private:
        int n_steps;
        int n_steps_terminate; //after this many steps, the game will be won by MIN.

        void set_vertex_type(int vertex, vertex_type type);
        void force_vertex_type(int vertex, vertex_type type);
        void force_edges(int vertex, int e1, int e2);
        void set_edges(int vertex, int e1, int e2);
        bool increment_min_strat(std::vector<bool> &s);
        bool increment_max_strat(std::vector<bool> &s);
        std::vector<bool> bruteforce_max(std::vector<bool> s);
};


std::ostream& operator<<(std::ostream& os, const std::vector<bool> &vec);


std::ostream& operator<<(std::ostream& os, const std::vector<double> &vec);

std::ostream& operator<<(std::ostream& os, const std::vector<int> &vec);