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

bool PathIsFeasible(Polygon_2 initialPolygon,Point_2 v3, Point_2 v1,int v3It,int v1It, int eSourceIt,int eTargetIt) {
    // cout << "moving " << v3 << "inside e " << endl; 
    // for (int i = 0; i < path.size(); i++) {
    //     // vertex_iterator = findIterator()
    // }
    // initialPolygon.erase(initialPolygon.begin() + v3It);
    // initialPolygon.insert(initialPolygon.begin() + eSourceIt , v3);
    // initialPolygon.erase(initialPolygon.begin() + v1It);
    // initialPolygon.insert(initialPolygon.begin() + eSourceIt, v1);
    // cout << "moving path to edge makes polygon :" << endl; 
    // printPolygonEdges(initialPolygon); 
    // if (initialPolygon.is_simple()) {
    //     return true;
    // }
    // return false; 
}