#pragma once

#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>

extern std::random_device rd;
extern std::default_random_engine rng;



/* returns a random permutation of numbers {start, start+1, start+2... end-1}
    note that generating ranges with large integers will put a huge vector into memory
*/
std::vector<int> permutation_in_range(int start, int end);