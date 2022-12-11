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


    vector<char*> argvsFirstAssign;
    char str0[] = "project-algo"; 
    char str1[] = "polygon"; 
    char str2[] = "-i"; 
    string str3 = argv[2]; 
    char str12[] = "r";
    char str4[] = "-o";
    char str5[] = "polygon"; 
    char str6[] = "-algorithm";
    char str7[] = "incremental"; 
    char str8[] = "-edge_selection";
    char str9[] = "3";            /// max area: 3, min area: 2 , rand: 1 
    char str10[] = "-initialization";
    char str11[] = "1a"; 
    char* argv1[] = {str1,str2,str12,str4,str5,str6,str7,str8,str9,str10,str11};  

    string algorithm = argv[6]; 

    if (algorithm == "local_search") {
        for (int i = 0; i < 11; i++) {
            argvsFirstAssign.push_back(argv1[i]); 
        }
    } else {
        for (int i = 0; i < 12; i++) {
            argvsFirstAssign.push_back(argv1[i]); 
        }
    }
    
    Polygon_2 initialPolygon;
    initialPolygon = test_polyg(argvsFirstAssign.size(),argvsFirstAssign, str3);   //getting the polygon from 1st assignment by calling test_polyg function  

    double initialPolygonArea = abs(initialPolygon.area());  
    Polygon_2 initialConvexHull = calc_convex_hull(initialPolygon); 
    double initialCHarea = abs(initialConvexHull.area());
    double initialRatio = initialPolygonArea/initialCHarea; 
  
    std::vector<pair<Point_2,int>> vertex_iterators;                //initializing a vector of pairs that holds each vertex of polygon and the 
    reform(vertex_iterators,initialPolygon);                         //correct iterator to acquire it.
 

    int L = stoi(argv[8]); 
    string area_polygonization = argv[9]; 
    bool max_area_polygonization,min_area_polygonization; 
    if (area_polygonization == "-max") {
        max_area_polygonization = true; 
        min_area_polygonization = false; 
    }else if (area_polygonization == "-min") {
        max_area_polygonization = false; 
        min_area_polygonization = true; 
    }
    
    int total_modifications = 0; 

    cout << "algorithm is " << algorithm << endl; 
    
    if (algorithm == "local_search") {
    clock_t start, end;
    start = clock();
    float DA;
    string fs(argv[11]); 
    float threshold = stof(fs);
    if (max_area_polygonization) {
        DA = 5.0; 
    }else {
        DA = 15.0;   
    }
    vector<pair<Segment_2,vector<Point_2>>> T; 
    vector<float>T_areas; 

    while(DA >= threshold) {   
        float initialArea = abs(initialPolygon.area()); 
        vector <Point_2> path;  
        bool has_solution = false; 
        int breakpoint = initialPolygon.vertices().size() - L ; 
        for (const Segment_2& e :initialPolygon.edges()) {            //for each edge of polygon 
            int eSourceIt = findIterator(vertex_iterators,e.source());   //finding iterator of u1 vertex of blue edge    
            for (const Point_2& vertex: initialPolygon.vertices()) {        
                int vertex_iterator = findIterator(vertex_iterators,vertex); 
                if (vertex_iterator > breakpoint) {                        
                    path.clear(); 
                    break; 
                }
                for (int i = 1; i <= L; i++) {
                    for (int j = 0 ; j < i; j++) {
                        path.push_back(vertex_iterators[j+vertex_iterator].first);      //creating the path to be inserted depending on L parameter 
                    }                                                                   // L increases its time so we get every possible path from L=1,L=2,L=3 etc.
                                                                                        //depending on the value of L. 
                reverse(path.begin(),path.end());                                       
    
                Polygon_2 testPolygon = initialPolygon;                 //using a testPolygon so initialPolygon will stay the same. if path solution is feasible
                                                                        //then must apply changes to initialPolygon as well 
                for (int i = 0; i < path.size(); i++) {                   ///inserting the path to edge  
                    vertex_iterator = findIterator(vertex_iterators,path[i]); 
                    testPolygon.erase(testPolygon.begin() + vertex_iterator);
                    testPolygon.insert(testPolygon.begin() + eSourceIt , path[i]);
                }

                float new_area = abs(testPolygon.area());   
                initialArea = abs(initialPolygon.area());                
                int remainsSimple = testPolygon.is_simple(); 
                bool pathIsFeasibleSol = feasibleSolution(new_area,initialArea,remainsSimple, T, T_areas,path,max_area_polygonization,
                                        min_area_polygonization,e); 
                
                testPolygon.clear();           //clearing the testPolygon for next iteration
                path.clear();                   //clearing the path for next iteration 
                }
        
            }
        }
        if (T_areas.size() != 0) {
            int it; 
            if (max_area_polygonization) {
                it = distance(T_areas.begin(),max_element(T_areas.begin(), T_areas.end()));
            }else if (min_area_polygonization) {
                it = distance(T_areas.begin(),min_element(T_areas.begin(), T_areas.end()));
            }

            vector <Point_2> winningPath; 
            for (int j = 0; j < T[it].second.size(); j++) {
                winningPath.push_back(T[it].second[j]); 
            }

            int eSourceIt = findIterator(vertex_iterators,T[it].first.source()); 
            Polygon_2 potentialPolygon = initialPolygon; 
            for (int i = 0; i < winningPath.size(); i++) {                   ///inserting path to edge 
                int vertex_iterator = findIterator(vertex_iterators,winningPath[i]); 
                potentialPolygon.erase(potentialPolygon.begin() + vertex_iterator);
                potentialPolygon.insert(potentialPolygon.begin() + eSourceIt , winningPath[i]);
            }
            
            if (max_area_polygonization) {
                initialArea = abs(initialPolygon.area());
                DA = abs(potentialPolygon.area()) - initialArea;  // recalculate DA
                initialPolygon = potentialPolygon;                 //new Polygon with the best solution applied to it 
                reform(vertex_iterators,initialPolygon);            //reform vertex iterators so we get the correct iterator for each vertex for the new Polygon
                total_modifications++; 
            }else if (min_area_polygonization) {
                initialArea = abs(initialPolygon.area());
                DA = abs(initialPolygon.area()) - abs(potentialPolygon.area());
                initialPolygon = potentialPolygon; 
                reform(vertex_iterators,initialPolygon); 
                total_modifications++; 
            }

        T.clear(); 
        T_areas.clear();  
       }else {
            break;                //if there are no possible solutions to be applied to polygon, break out of the loop and end the proccess 
       }
        
    }
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    time_taken *=1000;


    for(const Point_2& vertex: initialPolygon.vertices()) {
        cout << vertex << endl; 
    }
    for(const Segment_2& edge: initialPolygon.edges()) {
        cout << edge << endl; 
    }
    cout << "Algorithm: " << algorithm << endl; 
    cout << "area_initial: " << initialPolygonArea << endl; 
    cout << "area: " << abs(initialPolygon.area()) << endl; 
    cout << "ratio_initial: "<< initialRatio << endl; 

    Polygon_2 finalCH = calc_convex_hull(initialPolygon); 
    double finalCHarea = abs(finalCH.area()); 
    cout << "ratio: " << abs(initialPolygon.area())/finalCHarea << endl; 
    cout << "Construction time: " << fixed << round(time_taken) << setprecision(5);
    cout << " ms " << endl;    



    }else if(algorithm =="simulated_annealing"){

        clock_t start, end;
        start = clock();

        string ann = argv[11];
        string m = argv[9];

        Polygon_2 ch;
        CGAL::convex_hull_2( initialPolygon.begin(), initialPolygon.end(), std::back_inserter(ch));
        int n=initialPolygon.size();

        Polygon_2 new_p = simulated_annealing(initialPolygon,ch,n,m,ann,L); 

        end = clock();

        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        time_taken *=1000;
        
        std::cout << "Optimal Area Polygonization"<< std::endl;
    
        for(auto it = new_p.begin(); it!= new_p.end();++it){
            std::cout << *it << std::endl;
        }

        for(const Segment_2& e: new_p.edges()){
            std::cout << e << std::endl;
        }

        Polygon_2 chn;
        CGAL::convex_hull_2( initialPolygon.begin(), initialPolygon.end(), std::back_inserter(chn));

        std::cout << "Algorithm: "<< algorithm << m << std::endl;
        std::cout << "area initial: " << initialPolygon.area() << std::endl;
        std::cout << "area: " << new_p.area() << std::endl;
        std::cout << "ratio_initial: " << initialPolygon.area()/ch.area() << std::endl;
        std::cout << "ratio: " << new_p.area()/chn.area() << std::endl;
        cout << "Construction time: " << fixed << round(time_taken) << setprecision(5);
        cout << " ms " << endl; 
    }
    
    return 0; 
}
