#include <CGAL/convex_hull_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <algorithm>   
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Origin.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Point_2.h>
#include <CGAL/Origin.h>
#include <CGAL/random_selection.h>
#include <CGAL/Polygon_set_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef std::vector<Point_2> Points;
typedef K::Intersect_2 Intersect_2;
typedef CGAL::Triangle_2<K> Triangle_2;
typedef CGAL::Random_points_in_square_2<Point_2> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;
typedef CGAL::Polygon_set_2<K> Polygon_set_2;


using std::cout; using std::endl;
using namespace std;

void getargv(std::vector<char*>&); 

Segment_2 findPriorEdge(Point_2,Polygon_2); 
Segment_2 findNextEdge(Point_2,Polygon_2); 
void printPolygonEdges(Polygon_2);

Polygon_2 simulated_annealing(Polygon_2 ,Polygon_2 , int , string , string , int);
bool is_intersected_P(Segment_2 , Polygon_2 );
bool is_intersected_S(Segment_2 , Segment_2 );
int position(std::vector<std::pair<Point_2,int>> , Point_2 );
void change_after(std::vector<std::pair<Point_2,int>>& ,Polygon_2 );
Point_2 find_maxx(Polygon_2 );
Point_2 find_minx(Polygon_2);
