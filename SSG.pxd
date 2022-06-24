cdef extern from "SSG.cpp":
    pass

# Declare the class with cdef
cdef extern from "include/SSG.h" namespace "game":
    cdef cppclass SSG:
        Rectangle() except +
        Rectangle(int, int, int, int) except +
        int x0, y0, x1, y1
        int getArea()
        void getSize(int* width, int* height)
        void move(int, int)

cdef cppclass SSG:
        int n;

        std::vector<std::vector<int>> outgoing_edge
        std::vector<std::vector<int>> incoming_edges
        std::vector<vertex_type> type
        std::vector<int> min_vertices
        std::vector<int> max_vertices
        std::vector<int> ave_vertices
        int min_sink_vertex = -1
        int max_sink_vertex = -1
        int token
        std::vector<bool> strat
        double tolerance
        const int c
        double beta
   
        SSG(int n_vertices) except +
        void set_vertex(int vertex, vertex_type type, int e1, int e2)


        void start(int starting_vertex, std::vector<bool> strategy)
        int step(int n_steps = 1)
        int play(int starting_vertex, std::vector<bool> strategy)
        double play_n(int starting_vertex, std::vector<bool> strategy, int n_trials)

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
        static SSG random_game_equal_split(int n);
        static SSG random_game_mod(int n);
        static std::vector<bool> random_strategy(int n);
        static SSG read_game_file(std::ifstream &file);
        static std::vector<bool> read_strategy_file(std::ifstream &file);
        static bool probs_match(const std::vector<double> &p1, const std::vector<double> &p2, double tolerance);


};