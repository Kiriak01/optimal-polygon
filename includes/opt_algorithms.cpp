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

void getargv(std::vector<char*>& argvsFirstAssign) {
    
    
}
