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
        {"local_search", "incremental"} , {"local_search", "convexhull"}
    };

    vector <std::string> algo_results;
    double max_ratio, min_ratio, cut_off;
    int points_amount; 
    
    vector <pair <std::string , int>> algo_points;
    vector <pair <double, double>> algo_res; 
    

    for (int i = 0; i < algorithms.size(); i++) {
        string opt_algorithm = algorithms[i].first;
        string init_algorithm = algorithms[i].second; 
        string total_algorithm = opt_algorithm + "-" + init_algorithm; 
        for (int j = 0; j < files.size(); j++) {
            points_amount = getSumPoints(files[j]); 
            cut_off  = points_amount * 500; //(ms)
            max_ratio = opt_local_search(files[j], opt_algorithm, 1, "-max", "0.1",
                                            init_algorithm,  2, "2b", "global", cut_off); 
            min_ratio = opt_local_search(files[j], opt_algorithm, 1, "-min", "0.1",
                                           init_algorithm,  2, "2b", "global", cut_off); 
        
            algo_points.push_back(make_pair(total_algorithm,points_amount));
            algo_res.push_back(make_pair(max_ratio,min_ratio)); 
            
        }
    }

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

    getBestResults(algo_points,algo_res,max_ratio_algo,min_ratio_algo); 

    printResultsBoard(max_ratio_algo,min_ratio_algo); 

    // for (int i = 0; i < max_ratio_algo.size(); i++) {
    //     pair <std::string, double> p;
    //     p = max_ratio_algo[i].first; 
    //     int points = max_ratio_algo[i].second;
    //     cout << p.first << " , " << p.second << " , " << points << endl; 
    // }
   

    return 0; 
}
