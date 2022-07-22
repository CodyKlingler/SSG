#pragma once

#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <set>

extern std::random_device rd;
extern std::default_random_engine rng;


/* returns a random permutation of numbers {start, start+1, start+2... end-1} */
std::vector<int> permutation_in_range(int start, int end);

int random_in_range_exclude(int start, int end, std::set<int> exclude);
