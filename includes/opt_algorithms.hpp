#include <CGAL/convex_hull_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <algorithm>   
#include <vector>
#include <cmath>
#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Origin.h>




typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef std::vector<Point_2> Points;
typedef Kernel::Intersect_2 Intersect_2;
typedef CGAL::Triangle_2<Kernel> Triangle_2;



typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef std::vector<Point_2> Points;
using std::cout; using std::endl;

Segment_2 findPriorEdge(Point_2,Polygon_2); 
Segment_2 findNextEdge(Point_2,Polygon_2); 
void printPolygonEdges(Polygon_2); 
bool feasibleSolution(float,float,int,std::vector<std::pair<Segment_2,std::vector<Point_2>>>&, std::vector<float>&, std::vector<Point_2>, bool, bool, Segment_2);


