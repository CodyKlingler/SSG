#include "Strategy.h"


Strategy::Strategy(){
    n=-1;
    strategy=NULL;
}

Strategy::Strategy(SSG game){
    int vertices = game.n;
    strategy = new int[vertices]{0};
}


Strategy::Strategy(int vertices){
    strategy = new int[vertices]{0};
}