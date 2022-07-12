#pragma once

#include "include/permute.h"
#include <random>


std::random_device rd = std::random_device{};
std::default_random_engine  rng = std::default_random_engine{rd()};

std::vector<int> n_vec(0);


//returns a random permutation of numbers {start, start+1, start+2... end-1}
std::vector<int> permutation_in_range(int start, int end){
    if(start > end || start < 0){
        std::vector<int> v;
        return v;
    }

    //extend size of n_vec if necessary
    while(n_vec.size()<end){
        n_vec.push_back((int)n_vec.size());
    }

    std::vector<int> ret(n_vec.begin()+start, n_vec.begin()+end);

    std::shuffle(std::begin(ret), std::end(ret), rng);

    return ret;
}