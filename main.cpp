#include <iostream>
#include <algorithm>   
#include <vector>
#include "project-algo/test_polyg.cpp"
#include "includes/opt_algorithms.hpp"
#include "project-algo/includes/algorithms.hpp"
#include "project-algo/includes/algorithms.cpp"
#include <string> 
#include <fstream>
using namespace std; 


int main(int argc, char* argv[]) {


    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "simulated_annealing", 3, "-max", "0.1", "incremental",  3, "1a"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "local_search", 1, "-max", "0.1", "incremental",  3, "1b", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "simulated_annealing", 3, "-max", "0.1", "convexhull",  1, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000040.instance", "local_search", 5, "-min", "0.1", "incremental",  2, "1a"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "simulated_annealing", 3, "-max", "0.1", "incremental",  3, "1a"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000050.instance", "local_search", 8, "-max", "0.1", "incremental",  3, "1a"); 
    
    opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000050.instance", "local_search", 2, "-max", "0.1", "convexhull",  2, "2b", "global"); 


    // for (int i = 1; i < 11; i++) {
    //     // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000030.instance", "local_search", i, "-max", "0.1", "convexhull",  3, "2b", "global"); 
    //     // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000010.instance", i, "-min", "0.1", "convexhull",  2, "1b"); 
    // }
    // paikse me tis parametrous. epishs ftiakse vectors gia max kai gia min antistoixa kai des tis megalyteres kai mikroteres times 
    // kai pote autes prokyptoyn(gia poies parametrous kai gia poia arxeia).  

    return 0; 
}
