#include "include/permute.h"

std::random_device rd = std::random_device{};
std::default_random_engine  rng = std::default_random_engine{rd()};

//returns a random permutation of numbers {start, start+1, start+2... end-1}
std::vector<int> permutation_in_range(int start, int end){
    if(start >= end || start < 0){
        std::vector<int> v;
        return v;
    }

    const int n_in_range = end-start;

    std::vector<int> ret(n_in_range);

    for(int i = 0; i<(n_in_range); i++){
        ret[i] = start + i;
    }

    std::shuffle(std::begin(ret), std::end(ret), rng);

    return ret;
}

//returns a random permutation of numbers {start, start+1, start+2... end-1}
std::set<int> range(int start, int end){
    if(start >= end){
    }

    const int n_in_range = end-start;

    std::set<int> ret;

    for(int i = 0; i<(n_in_range); i++){
        ret.insert(start+i);
    }

    return ret;
}

//return random integer in the set {start, start+1, start+2... end-1} - exclude
int random_in_range_exclude(int start, int end, std::set<int> exclude){
    if(end <= start)
        return start;

    int n;

    auto set = range(start, end);

    for(auto it = exclude.rbegin(); it!=exclude.rend(); it++){
        set.erase(*it);
    }

    if(set.size()){

        int index = random()%set.size();
        auto it = set.begin();
        for(int i = 0; i<index; i++){
            it++;
        }
        return *it;
    }
    

    return 0;
}