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


    vector< pair <std::string, std::string>> algorithms = {
        {"local_search", "incremental"} , {"local_search", "convexhull"},
        {"simulated_annealing", "incremental"} //, {"simulated_annealing", "convexhull"}  sim an with ch too slow for now uncomment later 
    };

    vector <std::string> annealing_step = { {"local"}, {"global"}, {"subdivisor"}};   //ta prostetoyme sto for otan valoyme cutoff ston annealing 

    vector <std::string> algo_results;
    double cut_off;
    int points_amount; 
    
    vector <pair <std::string , int>> algo_points;
    vector <pair <double, double>> algo_res; 
    string opt_algorithm, init_algorithm, total_algorithm; 
    
    
    // vector <pair <pair <std::string , double > , int >> max_bound; 
    // vector <pair <pair <std::string , double > , int >> min_bound; 
    vector<double> test_max_bound; 
    pair <std::string, double> testPmax, pMax, pMin;
    int points; 
    string algo;
    double score_max; 
    
    double max_ratio, min_ratio ,max_bound, min_bound; 
    

    for (int i = 0; i < algorithms.size(); i++) {
        opt_algorithm = algorithms[i].first;
        init_algorithm = algorithms[i].second; 
        total_algorithm = opt_algorithm + "-" + init_algorithm; 
        int L; 
        if (opt_algorithm == "local_search"){
            L = 1;
        }else {
            L = 50;  
        }
        max_bound = 1.0; 
        min_bound = 0.0; 
        for (int j = 0; j < files.size(); j++) {
            points_amount = getSumPoints(files[j]); 
            cut_off  = points_amount * 500; //(ms)
            max_ratio = optimize_polygon(files[j], opt_algorithm, L, "-max", "0.1",
                                            init_algorithm,  2, "2b", "local", cut_off); 
            min_ratio = optimize_polygon(files[j], opt_algorithm, L, "-min", "0.1",
                                           init_algorithm,  2, "2b", "local", cut_off); 
        
            algo_points.push_back(make_pair(total_algorithm,points_amount));
            algo_res.push_back(make_pair(max_ratio,min_ratio));
        }
    }

    vector <pair <pair <std::string, double>, int>> MAX_bounds; 
    vector <pair <pair <std::string, double>, int>> MIN_bounds; 
    createBoundVectors(algo_points,algo_res,MAX_bounds,MIN_bounds);  


    vector <pair <pair <std::string , double > , int >> max_ratio_algo; 
    vector <pair <pair <std::string , double > , int >> min_ratio_algo; 
    for (int i = 0; i < algorithms.size(); i++) {
        for (int j = 0; j < files.size(); j++) {
            points_amount = getSumPoints(files[j]); 
            pair <std::string , double> p;
            p.first = algorithms[i].first + "-" + algorithms[i].second;
            p.second = 0.0;
            max_ratio_algo.push_back(make_pair(p,points_amount));
            min_ratio_algo.push_back(make_pair(p,points_amount));
        }
    }



    getBestResults(algo_points,algo_res,max_ratio_algo,min_ratio_algo,files); 

    printResultsBoard(max_ratio_algo,min_ratio_algo); 
    

    return 0; 
}
