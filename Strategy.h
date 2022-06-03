#pragma once
#include "SSG.h"

class Strategy{
    public:

        int n;

        int* strategy;
        Strategy();
        Strategy(int n_vertices);
        Strategy(const SSG game);
        Strategy(const SSG game, const Strategy min_strategy, const Strategy max_strategy);

        inline int operator [] (int i){return strategy[i];}
};