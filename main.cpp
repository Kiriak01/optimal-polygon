#include <iostream>
#include <algorithm>   
#include <vector>
#include "includes/opt_algorithms.hpp"
#include "project-algo/test_polyg.cpp"
#include "project-algo/includes/algorithms.hpp"
#include "project-algo/includes/algorithms.cpp"
#include <string> 
#include <fstream>
using namespace std; 

int main(int argc, char* argv[]) {

    //run ./optimal_polygon -i euro-night-0000010.instance -o optimal_polygon -algorithm local_search -L 2 -max -threshold 0.1 
    // to push git push -f optimal-project-algo master

    vector<char*> argvsFirstAssign;
    char str0[] = "project-algo"; 
    char str1[] = "polygon"; 
    char str2[] = "-i"; 
    char str3[] = "euro-night-0000010.instance"; 
    char str4[] = "-o";
    char str5[] = "polygon"; 
    char str6[] = "-algorithm";
    char str7[] = "incremental"; 
    char str8[] = "-edge_selection";
    char str9[] = "3";            /// max area: 3, min area: 2 , rand: 1 
    char str10[] = "-initialization";
    char str11[] = "1a"; 
    char* argv1[] = {str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11};  
    for (int i = 0; i < 11; i++) {
        std::cout << argv1[i] << "\n";
        argvsFirstAssign.push_back(argv1[i]); 
    }    
    // for (int i =0; i < argvsFirstAssign.size(); i ++) {
    //     cout << "thesi: " << i << " , " << argvsFirstAssign[i] << endl; 
    // }
    
    Polygon_2 initialPolygon;
    initialPolygon = test_polyg(argvsFirstAssign.size(),argvsFirstAssign); 
    double initialPolygonArea = abs(initialPolygon.area());  
    cout << "initial poly has area " << initialPolygonArea << endl; 
    
    cout << "is simple: " << initialPolygon.is_simple() << endl;


    std::vector<pair<Point_2,int>> vertex_iterators; 
    reform(vertex_iterators,initialPolygon); 
    for(int i = 0; i < vertex_iterators.size(); i++)
    {
        cout << vertex_iterators[i].first << "," << vertex_iterators[i].second << endl; 
    }

    string algorithm = argv[6]; 
    int L = stoi(argv[8]); 
    string area_polygonization = argv[9]; 
    string fs(argv[11]); 
    float threshold = stof(fs);
    cout << "algorithm is : " << algorithm << endl; 
    cout << "L is: " << L << endl; 
    cout << "area polygonization is: " << area_polygonization << endl; 
    cout << "threshold is: " << threshold << endl; 

    if (algorithm == "local_search") {
        cout << "local search algo " << endl; 

        // float DA = 1.0; 
        // while(DA >= threshold) {
        //     for(const Segment_2& e: initialPolygon.edges()){
        //         for (int i = 0; i < L; i++) {
        //             ///some code here for the local search algorithm 
        //         }
                
        //     }
        // }
        //
        //


    }


    
    
    return 0; 
}
