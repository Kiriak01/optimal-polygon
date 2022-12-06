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
    // initialPolygon = test_polyg(argvsFirstAssign.size(),argvsFirstAssign);   //uncomment to execute 1st assignment code and get initial polygon 
    initialPolygon.push_back(Point_2(0,0)); 
    initialPolygon.push_back(Point_2(4,4)); 
    initialPolygon.push_back(Point_2(6,6)); 
    initialPolygon.push_back(Point_2(8,6)); 
    initialPolygon.push_back(Point_2(7,4)); 
    initialPolygon.push_back(Point_2(10,2)); 
    initialPolygon.push_back(Point_2(7,0)); 

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
        printPolygonEdges(initialPolygon); 
        // for (int i = 0; i < vertex_iterators.size(); i++) {
        //     cout << vertex_iterators[i].first << " , " << vertex_iterators[i].second << endl; 
        // }
        // Segment_2 e(Point_2(7,0),Point_2(0,0)); 
        // int eSourceIt, eTargetIt = -1; 
        // eSourceIt = findIterator(vertex_iterators,e.source());  
        // eTargetIt = findIterator(vertex_iterators,e.target()); 

        // cout << "for edge " << e << endl; 
        // vector <Point_2> myPath; 
        // myPath.push_back(Point_2(6,6));
        // // myPath.push_back(Point_2(8,6));
        // // myPath.push_back(Point_2(7,4));
        // reverse(myPath.begin(),myPath.end()); 
        // cout << "path is: " << endl; 
        // for (int i = 0; i < myPath.size(); i ++) {
        //     cout << myPath[i] << endl; 
        // }
        

    float DA = 5.0; 

    while (DA >= threshold) {       //DA
        // cout << "DA , THRESHOLD " << DA << "," << threshold << endl; 
        vector<pair<Segment_2,vector<Point_2>>> T; 
        vector<float>T_areas; 
        float initialArea = abs(initialPolygon.area()); 
        vector <Point_2> path;  
        bool has_solution = false; 
        int breakpoint = initialPolygon.vertices().size() - L ; 
        for (const Segment_2& e :initialPolygon.edges()) {
            int eSourceIt = findIterator(vertex_iterators,e.source());   //finding iterator of u1 vertex of blue edge    
            for (const Point_2& vertex: initialPolygon.vertices()) {           
                int vertex_iterator = findIterator(vertex_iterators,vertex); 
                if (vertex_iterator > breakpoint) {
                    path.clear(); 
                    break; 
                }
                for (int i = 1; i <= L; i++) {
                    for (int j = 0 ; j < i; j++) {
                        path.push_back(vertex_iterators[j+vertex_iterator].first); 
                    }
                reverse(path.begin(),path.end()); 
    
                Polygon_2 testPolygon = initialPolygon;                 //using a testPolygon so initialPolygon will stay the same. if path solution is feasible
                                                                        //then must apply changes to initialPolygon as well 
                for (int i = 0; i < path.size(); i++) {                   ///inserting path to edge 
                    vertex_iterator = findIterator(vertex_iterators,path[i]); 
                    testPolygon.erase(testPolygon.begin() + vertex_iterator);
                    testPolygon.insert(testPolygon.begin() + eSourceIt , path[i]);
                }

                float new_area = abs(testPolygon.area());                  
                int remainsSimple = testPolygon.is_simple(); 
                if ((new_area > initialArea) && (remainsSimple)) {         //checking if transfer is feasible
                        cout << "moving path to " << e << " increases area and retains simplicity " << endl; 
                        cout << "new area , initialArea : " << new_area << ", " << initialArea << endl; 
                        T.push_back(make_pair(e , path));         
                        T_areas.push_back(new_area); 
                        has_solution = true; 
                    }else {
                        cout << "not feasible solution because areas: " << new_area << "," << initialArea << " and also simple: "<< remainsSimple << endl; 
                    }
                testPolygon.clear();
                testPolygon = initialPolygon;
                path.clear(); 
                }
        
            }
        }
        if (has_solution) {
            cout << "MPAINWWWWWWWWWWWWWWWWWWWWWWWWWWWW" << endl; 
            cout << "printing T solutions" << endl; 
            for (int i =0; i<T.size(); i++) {
                cout << "Area: " << T_areas[i] << " from edge " << T[i].first << " with path: "; 
                for (int j = 0; j < T[i].second.size(); j++) {
                    cout << T[i].second[j] << endl; 
                }
            }
            cout << "polygon back to its initial state" << endl; 
            printPolygonEdges(initialPolygon); 
            int it = distance(T_areas.begin(),max_element(T_areas.begin(), T_areas.end()));
            cout << "max area is " << T_areas[it] << "with iterator " << it << endl; 

            cout <<"max edge: "<< T[it].first  << " with path " << endl; 
            vector <Point_2> winningPath; 
            for (int j = 0; j < T[it].second.size(); j++) {
                winningPath.push_back(T[it].second[j]); 
            }


            int eSourceIt = findIterator(vertex_iterators,T[it].first.source()); 
            cout << "source edge is " << T[it].first.source() << " ,  eSOurce it is " << eSourceIt << endl; 
            Polygon_2 potentialPolygon = initialPolygon; 
            

            for (int i = 0; i < winningPath.size(); i++) {                   ///inserting path to edge 
                int vertex_iterator = findIterator(vertex_iterators,winningPath[i]); 
                potentialPolygon.erase(potentialPolygon.begin() + vertex_iterator);
                potentialPolygon.insert(potentialPolygon.begin() + eSourceIt , winningPath[i]);
            }
            cout << "potential polygon " << endl; 
            printPolygonEdges(potentialPolygon); 
            cout << "initial polygon " << endl;
            printPolygonEdges(initialPolygon); 
            DA = abs(potentialPolygon.area()) - abs(initialPolygon.area());  
            cout << "DA IS " << DA << endl; 

            if (T_areas[it] > initialArea) {
                initialArea = T_areas[it]; 
                initialPolygon = potentialPolygon; 
                cout << "intialPolyon = potentialPolyg" << endl; 
                printPolygonEdges(initialPolygon) ;
                reform(vertex_iterators,initialPolygon); 
            }
        }

    }



    }
    
    return 0; 
}
