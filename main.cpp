#include <iostream>
#include <algorithm>   
#include <vector>
#include "project-algo/test_polyg.cpp"
#include "includes/opt_algorithms.hpp"
#include "project-algo/includes/algorithms.hpp"
#include "project-algo/includes/algorithms.cpp"
#include <string> 
#include <fstream>
#include <dirent.h> 

using namespace std; 


int main(int argc, char* argv[]) {


    vector<std::string> files; 

    DIR *dir;
    struct dirent *ent;
    string backslash = "/";
    if ((dir = opendir (argv[2])) != NULL) {            //getting all the set point files from directory and storing them into a vector
    while ((ent = readdir (dir)) != NULL) {
        char * file = strdup(ent->d_name); 
        if (strcmp(file,".") == 0 || strcmp(file,"..") == 0) {
            continue;
        } else {
            string path = argv[2]; 
            string filename = file; 
            string full_path = path + backslash + filename; 
            files.push_back(full_path); 
        }
    }
    closedir (dir);
    } else {
    perror ("");              //error at opening directory 
    return EXIT_FAILURE;
    }

    // cout << "files " << endl; 
    for (int i = 0; i < files.size(); i++) {
        int points_amount = getSumPoints(files[i]); 
        // cout << files[i] << endl; 
    
        // opt_local_search(files[i], "local_search", 3, "-max", "0.1", "incremental",  3, "1b", "global"); 
    }

    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "simulated_annealing", 3, "-max", "0.1", "incremental",  3, "1a"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "local_search", 1, "-max", "0.1", "incremental",  3, "1b", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "simulated_annealing", 3, "-max", "0.1", "convexhull",  1, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000040.instance", "local_search", 5, "-min", "0.1", "incremental",  2, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "simulated_annealing", 3, "-max", "0.1", "incremental",  3, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000050.instance", "local_search", 8, "-max", "0.1", "incremental",  3, "1a", "global"); 
    
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000010.instance", "local_search", 5, "-max", "0.1", "convexhull",  2, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "local_search", 5, "-max", "0.1", "convexhull",  2, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000030.instance", "local_search", 5, "-max", "0.1", "convexhull",  2, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000040.instance", "local_search", 5, "-max", "0.1", "convexhull",  2, "1a", "global"); 
    
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000060.instance", "local_search", 2, "-max", "0.1", "incremental",  2, "1b", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000070.instance", "local_search", 2, "-max", "0.1", "incremental",  2, "1b", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000080.instance", "local_search", 2, "-max", "0.1", "incremental",  2, "1b", "global"); 
    
    
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000010.instance", "local_search", 5, "-min", "0.1", "incremental",  3, "1a", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000020.instance", "local_search", 5, "-min", "0.1", "incremental",  3, "1b", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000030.instance", "local_search", 5, "-min", "0.1", "incremental",  3, "2a", "global"); 
    
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000040.instance", "local_search", 2, "-min", "0.1", "incremental",  2, "2b", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000050.instance", "local_search", 2, "-min", "0.1", "incremental",  2, "2b", "global"); 
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000060.instance", "local_search", 2, "-min", "0.1", "incremental",  2, "2b", "global"); 
    opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000500.instance", "local_search", 2, "-min", "0.1", "incremental",  2, "2b", "global"); 


    //gia min = 100,200 shmeia ta kalytera apotelesmata dinei o local search me L=2, kai incremental me max polygonization(3 = edge selection)

    
    
    // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000100.instance", "local_search", 2, "-max", "0.1", "incremental",  2, "1b", "global"); 

    // for (int i = 1; i < 11; i++) {
    //     // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000030.instance", "local_search", i, "-max", "0.1", "convexhull",  3, "2b", "global"); 
    //     // opt_local_search("/home/george/Desktop/instances/data/images/euro-night-0000010.instance", i, "-min", "0.1", "convexhull",  2, "1b"); 
    // }
    // paikse me tis parametrous. epishs ftiakse vectors gia max kai gia min antistoixa kai des tis megalyteres kai mikroteres times 
    // kai pote autes prokyptoyn(gia poies parametrous kai gia poia arxeia).  

    return 0; 
}
