#pragma once

#include <iostream>
#include <vector>

#include "SSG.h"

SSG condon_game();

void test_condon_example();

void test_hoffman(int n_tests, int n_strats_per_game, int n_vertices);

void test_randomized_hoffman(int n_tests, int n_strats_per_game, int n_vertices);

void benchmark_SSG(int n_games, int n_vertices);
void benchmark_SSG(std::vector<SSG> games);

void benchmark_SSG_distribution(std::vector<SSG> games);
void benchmark_SSG_distribution2(std::vector<SSG> games, double min, double max, double increment);

bool test_correctness(int n_games, int n_strats_per_game, int n_vertices);

double test_stopping_constant(int n_games, int n_strats_per_game, int n_vertices);