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


using namespace std;
 

// export a game to a file
// export strategies to a file
// load game and strategies

//convert SSG::outgoing_edge 2 vector.


/*

std::chrono::time_point<std::chrono::system_clock> start, end;
  
    start = std::chrono::system_clock::now();
    std::cout << "f(42) = " << fibonacci(42) << '\n';
    end = std::chrono::system_clock::now();
  
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

*/



int main(){
    srand(time(NULL));
    std::cout << std::fixed << std::setprecision(3);

    //test_condon_example();  return 0;

    //test_hoffman_random(5, 100); return 0;


    /*
    ifstream myfile;
    myfile.open("this.txt");

    SSG g = SSG::read_game_file(myfile);

    vector<bool> s = {0,1,1,0,1,1};
    
    auto r = g.hoffman_karp(s);

    cout << r;

    return 0;
    */



   //TODO.. clean up debugging nonsense and make some benchmarking functions.
   //debugging stuff in hoffman() and probabilities()
    // -O3 optimization is on btw.
    // can speed hoffman up a lot by not iterating back through the vertices n^2 times
    


    for(int i = 5; i<10000; i++){
        auto start = std::chrono::system_clock::now();
        test_hoffman_random(3,i);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        cout << i << "\t" << elapsed_seconds.count() << endl;
    }

    return 0;
}