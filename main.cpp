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
    initialPolygon.push_back(Point_2(6,2)); 
    initialPolygon.push_back(Point_2(8,0)); 
    initialPolygon.push_back(Point_2(5,1)); 
    // initialPolygon.push_back(Point_2(8,4)); 
    // initialPolygon.push_back(Point_2(8,0)); 

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
        float initialArea = abs(initialPolygon.area()); //thelei me trigwnakia vlepe delos 

        printPolygonEdges(initialPolygon); 

        std::vector<pair<Segment_2,Point_2>> T; 
        std::vector<float>T_areas; 
        // Point_2 vertex(5,1);
        //while DA<=threshold {}
        for(const Segment_2& e :initialPolygon.edges()) {    
            cout << "checking for edge: " << e << endl; 
            int eSourceIt, eTargetIt = -1; 
            eSourceIt = findIterator(vertex_iterators,e.source());  
            eTargetIt = findIterator(vertex_iterators,e.target()); 
            cout << "segment iterators " <<  eSourceIt << " , " << eTargetIt << endl; 
            for(const Point_2& vertex: initialPolygon.vertices()) { 
            int vertex_iterator = -1;
            vertex_iterator = findIterator(vertex_iterators,vertex);
            Point_2 myVertex = vertex; 
            int myit = vertex_iterator; 
            cout << "adding vertex: " << myVertex << " to edge: " << e << endl; 
            if ((eSourceIt!=-1) && (eTargetIt!=-1)) {
                int maxIt = max(eSourceIt,eTargetIt);
                initialPolygon.erase(initialPolygon.begin() + vertex_iterator); 
                initialPolygon.insert(initialPolygon.begin() + maxIt , myVertex); 
                float new_area = abs(initialPolygon.area());                  //prepei na ginei me trigwnakia blepe delos 
                int remainsSimple = initialPolygon.is_simple(); 
                if ((new_area > initialArea) && (remainsSimple)) {
                    cout << "moving " << vertex << " to " << e << " increases area and retains simplicity " << endl; 
                    cout << "new area , initialArea : " << new_area << ", " << initialArea << endl; 
                    T.push_back(make_pair(e , myVertex));
                    T_areas.push_back(new_area); 
                }else {
                    cout << "not feasible solution because areas: " << new_area << "," << initialArea << " and also simple: "<< remainsSimple << endl; 
                }
                initialPolygon.erase(initialPolygon.begin()+maxIt);
                initialPolygon.insert(initialPolygon.begin()+myit,myVertex); 
                reform(vertex_iterators,initialPolygon);
            }else if (eSourceIt != -1) {
                cout << "MPHKA STO ESORUCE IT MONO TOY PWS EGINE AUTO REE" <<endl; 
                initialPolygon.erase(initialPolygon.begin() + vertex_iterator); 
                initialPolygon.insert(initialPolygon.begin() + eSourceIt , myVertex); 
                float new_area = abs(initialPolygon.area());                  //prepei na ginei me trigwnakia blepe delos 
                int remainsSimple = initialPolygon.is_simple(); 
                if ((new_area > initialArea) && (remainsSimple)) {
                    cout << "moving " << vertex << "increases area and retains simplicity " << endl; 
                    cout << "new area , initialArea : " << new_area << ", " << initialArea << endl; 
                    T.push_back(make_pair(e , myVertex));
                    T_areas.push_back(new_area); 
                }else {
                    cout << "not feasible solution because areas: " << new_area << "," << initialArea << " and also simple: "<< remainsSimple << endl; 
                }
                initialPolygon.erase(initialPolygon.begin()+eSourceIt);
                initialPolygon.insert(initialPolygon.begin()+vertex_iterator,vertex); 
                reform(vertex_iterators,initialPolygon);
            }
        }  
    }
            cout << "printing T solutions" << endl; 
            for (int i =0; i<T.size(); i++) {
                cout << T[i].first << "," << T[i].second << endl; 
                cout << "with area "<< T_areas[i] << endl; 
            }
            cout << "polygon back to its initial state" << endl; 
            printPolygonEdges(initialPolygon); 
            int it = distance(T_areas.begin(),max_element(T_areas.begin(), T_areas.end()));
            cout << T[it].first << "," << T[it].second << " with area: " << T_areas[it] << endl; 

            //prosthese to sto for gia na ginetai gia kathe akmh kai gia kathe vertex. molis ta katafereis kai doyleuei swsta vres pws tha ginei me path
            //epishs ti ginetai me to prwto vertex poy exei mono next edge?? . vlepw kai delos gia ypologismo emvadoy me trigwnakia 


    }
    
    return 0; 
}
