#include "opt_algorithms.hpp"
#include <iostream>
#include <vector>
#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/algorithm.h>

 
using namespace std;
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef std::vector<Point_2> Points;
typedef std::vector <Segment_2> Segments;
using std::cout; using std::endl;


Segment_2 findPriorEdge(Point_2 point, Polygon_2 polygon) {
    for (const Segment_2& edge: polygon.edges()) {
        if (edge.target() == point) {
            return edge; 
        }
    }
}
Segment_2 findNextEdge(Point_2 point, Polygon_2 polygon) {
    for (const Segment_2& edge: polygon.edges()) {
        if (edge.source() == point) {
            return edge; 
        }
    }
}

void printPolygonEdges(Polygon_2 polygon){
    for (const Segment_2& edge: polygon.edges()) {
        cout << edge << endl; 
    }
}

bool feasibleSolution(float new_area, float initialArea, int remainsSimple ,vector<pair<Segment_2
                    ,vector<Point_2>>>& T,vector<float>& T_areas, vector<Point_2> path, bool max_area_polygonization,
                    bool min_area_polygonization, Segment_2 e) 
{
    if (max_area_polygonization) { 
        if ((new_area > initialArea) && (remainsSimple)) {         //checking if transfer of path to e increases area and retains simplicity of polygon 
            // cout << "moving path to ncreases area and retains simplicity " << endl;         //if it does, then adding it to solutions vector T, and adding its area to vector T_areas
            // cout << "new area , initialArea : " << new_area << ", " << initialArea << endl; 
            T.push_back(make_pair(e , path));                          
            T_areas.push_back(new_area); 
            return true; 
            }else {
                return false; 
            }
    } else if (min_area_polygonization) {
        if ((new_area < initialArea) && (remainsSimple)) {             //checking if transfer is path to e decreases area and retains simplicity of polygon 
            // cout << "moving path to ncreases area and retains simplicity " << endl; 
            // cout << "new area , initialArea : " << new_area << ", " << initialArea << endl; 
            T.push_back(make_pair(e , path));         
            T_areas.push_back(new_area); 
            return true; 
            }else {
                return false; 
            }
    }
    return false; 
}
