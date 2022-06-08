
/*
#include "include/Strategy.h"
#include "include/SSG.h"


#define STRATEGY_print


Strategy::Strategy(){
    n=-1;
    strategy=NULL;
}


Strategy::Strategy(int vertices, bool* strat_arr){
    n = vertices;
    strategy = strat_arr;
}

Strategy::Strategy(SSG game){
    int vertices = game.n;
    strategy = new bool[vertices]{0};
}


Strategy::Strategy(int vertices){
    strategy = new bool[vertices]{0};
}


Strategy::~Strategy(){
    if(strategy){
        //delete[] strategy;
    }
}
*/