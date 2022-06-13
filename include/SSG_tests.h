#pragma once

#include <iostream>
#include <vector>

#include "SSG.h"

SSG condon_game();

void test_condon_example();

void test_hoffman(int n_tests, int n_strats_per_game, int n_vertices);

void test_randomized_hoffman(int n_tests, int n_strats_per_game, int n_vertices);

void benchmark_SSG(int n_games, int n_strats_per_game, int n_vertices);
