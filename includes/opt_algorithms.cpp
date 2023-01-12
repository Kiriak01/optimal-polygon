#include "opt_algorithms.hpp"
#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/algorithm.h>
#include <cstdlib>
#include <string>
#include <algorithm>

 
using namespace std;
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef std::vector<Point_2> Points;
typedef std::vector <Segment_2> Segments;
using std::cout; using std::endl;

//local search 
Segment_2 findPriorEdge(Point_2 point, Polygon_2 polygon) {
    for (const Segment_2& edge: polygon.edges()) {
        if (edge.target() == point) {
            return edge; 
        }
    }
    return Segment_2(Point_2(0,0), Point_2(0,0));    //NULL segment 
}
Segment_2 findNextEdge(Point_2 point, Polygon_2 polygon) {
    for (const Segment_2& edge: polygon.edges()) {
        if (edge.source() == point) {
            return edge; 
        }
    } 
    return Segment_2(Point_2(0,0), Point_2(0,0));    //NULL segment
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
                                                                    //if it does, then adding it to solutions vector T, and adding its area to vector T_areas  
            T.push_back(make_pair(e , path));                          
            T_areas.push_back(new_area); 
            return true; 
            }else {
                return false; 
            }
    } else if (min_area_polygonization) {
        if ((new_area < initialArea) && (remainsSimple)) {             //checking if transfer is path to e decreases area and retains simplicity of polygon 
            T.push_back(make_pair(e , path));         
            T_areas.push_back(new_area); 
            return true; 
            }else {
                return false; 
            }
    }
    return false; 
}

int getSumPoints(string file) {
    int pos = file.find(".");
    string sub = file.substr(0 , pos);
    string substr = sub.substr(sub.length() - 7); 
    int number = stoi(substr);
    return number ; 
}

Polygon_2 getInitialPolygon(string init_algo, int init_edge_selection, string algorithm, string filename) {
    vector<char*> argvsFirstAssign;
    char str0[] = "project-algo"; 
    char str1[] = "polygon"; 
    char str2[] = "-i"; 
    
    char str12[] = "r";
    char str4[] = "-o";
    char str5[] = "polygon"; 
    char str6[] = "-algorithm";
    
    char str7[50]; 

    if (init_algo == "incremental") {
        strncpy(str7, "incremental", 15);      
    }else if (init_algo == "convexhull"){
        strncpy(str7, "convexhull", 15);      
    }
    char str8[] = "-edge_selection";
    char str9[2];                   // max area: 3, min area: 2 , rand: 1 
    if (init_edge_selection == 2) {
        strncpy(str9, "2", 2);      //getting min polygon from 1st assignment
    }else {
        strncpy(str9, "3", 2);      //getting max polygon from 1st assignment
    }
    char str10[] = "-initialization";
    char str11[] = "1a"; 
    char* argv1[] = {str1,str2,str12,str4,str5,str6,str7,str8,str9,str10,str11};  

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
    initialPolygon = test_polyg(argvsFirstAssign.size(),argvsFirstAssign, filename);   //getting the polygon from 1st assignment by calling test_polyg function  
    return initialPolygon;
}


double optimize_polygon(string filename, string algorithm, int L, string area_polygonization, string threshold,
                      string init_algo, int init_edge_selection, string init_vertex_sort, string annealing, double cut_off) { 
    
    
    cout << "optimizing file: " << filename << endl;  
    Polygon_2 initialPolygon = getInitialPolygon(init_algo, init_edge_selection ,algorithm, filename);

    bool max_area_polygonization,min_area_polygonization; 
    string m; 
   

    if (algorithm == "local_search") {
        clock_t start, end;
        start = clock();

        double initialPolygonArea = abs(initialPolygon.area());  
        Polygon_2 initialConvexHull = calc_convex_hull(initialPolygon); 
        double initialCHarea = abs(initialConvexHull.area());
        double initialRatio = initialPolygonArea/initialCHarea; 
    
        std::vector<pair<Point_2,int>> vertex_iterators;                //initializing a vector of pairs that holds each vertex of polygon and the 
        reform(vertex_iterators,initialPolygon);                         //correct iterator to acquire it.
    
        if (area_polygonization == "-max") {
            max_area_polygonization = true; 
            min_area_polygonization = false; 
            m = "max"; 
        }else if (area_polygonization == "-min") {
            max_area_polygonization = false; 
            min_area_polygonization = true; 
            m = "min"; 
        }
        float DA;
        string fs = "0.1"; 
        float threshold = stof(fs);
        if (max_area_polygonization) {
            DA = 5.0; 
        }else {
            DA = 15.0;   
        }
        vector<pair<Segment_2,vector<Point_2>>> T; 
        vector<float>T_areas; 
        double total_cutoff = 0.0; 

        while(DA >= threshold) {   
            clock_t func_start, func_end;
            func_start = clock();
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
                    
                    if (testPolygon.size() != 0 && path.size() != 0) {
                        testPolygon.clear();           //clearing the testPolygon for next iteration
                        path.clear();                   //clearing the path for next iteration 
                    }
                }
            
                }
                func_end = clock();
                double func_time_taken = double(func_end - func_start) / double(CLOCKS_PER_SEC);
                func_time_taken *=15;
                total_cutoff+=func_time_taken;
                if (total_cutoff >= cut_off) {                          //if optimazition exceeds cutoff time, terminate execution 
                    if (max_area_polygonization){
                        return 0.0;
                    }else {
                        return 1.0; 
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

                initialArea = abs(initialPolygon.area());
                DA = abs(potentialPolygon.area()) - initialArea;  // recalculate DA
                initialPolygon = potentialPolygon;                 //new Polygon with the best solution applied to it 
                reform(vertex_iterators,initialPolygon);            //reform vertex iterators so we get the correct iterator for each vertex for the new Polygon
                
                if (T.size()!=0 && T_areas.size()!=0) {
                    T.clear(); 
                    T_areas.clear();  
                }
        }else {
                break;                //if there are no possible solutions to be applied to polygon, break out of the loop and end the proccess 
        }
            
        }
        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        time_taken *=1000;

        Polygon_2 finalCH = calc_convex_hull(initialPolygon); 
        double finalCHarea = abs(finalCH.area()); 
        double ratio = abs(initialPolygon.area()) / finalCHarea; 
        return ratio;

    } else {

        clock_t start, end;
        start = clock();
        Polygon_2 ch;
        
        CGAL::convex_hull_2( initialPolygon.begin(), initialPolygon.end(), std::back_inserter(ch));
        
        int n=initialPolygon.size();
        string m = area_polygonization;
        string ann = annealing; 
        bool exceed_cutoff = false; 
        Polygon_2 new_p = simulated_annealing(initialPolygon,ch,n,m,ann,L,cut_off,exceed_cutoff);  
        end = clock();

        
        
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        time_taken *=1000;

        Polygon_2 chn;
        CGAL::convex_hull_2( new_p.begin(), new_p.end(), std::back_inserter(chn));
    
        double ratio; 
        if (exceed_cutoff == true) {
            if (m=="-max") {
                ratio = 0.0;
                return ratio; 
            }else {
                ratio = 1.0; 
                return ratio; 
            }
        }else {
            ratio = abs(new_p.area())/abs(chn.area());
            return ratio; 
        }

        return ratio; 

    }
}; 


void getBestResults(std::vector <std::pair <std::string , int>> & algo_points,  std::vector <std::pair <double, double>> & algo_res , 
                    std::vector <std::pair <std::pair <std::string , double > , int >> & max_ratio_algo,
                    std::vector <std::pair <std::pair <std::string , double > , int >> & min_ratio_algo, std::vector<std::string> files) 
{
    string algorithm, algo;
    int points,new_points;
    double max_ratio, min_ratio , score_max, score_min;
    pair <std::string, double> pMax, pMin, new_p_max, new_p_min;

    for (int i = 0; i < algo_points.size(); i++) {
        algorithm = algo_points[i].first;
        points = algo_points[i].second; 
        max_ratio = algo_res[i].first;
        min_ratio = algo_res[i].second; 
        for (int j = 0; j < max_ratio_algo.size(); j++) {
            pMax = max_ratio_algo[j].first; 
            pMin = min_ratio_algo[j].first; 
            algo = pMax.first;
            score_max = pMax.second;
            score_min = pMin.second; 
            new_points = max_ratio_algo[j].second; 
            if ((algorithm == algo) && (points==new_points)) {
                pMax.second += max_ratio;
                pMin.second += min_ratio;
                new_p_max.first = algorithm;
                new_p_max.second = pMax.second;
                new_p_min.first = algorithm;
                new_p_min.second = pMin.second;
                max_ratio_algo[j].first = new_p_max;
                min_ratio_algo[j].first = new_p_min;
            }
        }
    }

    cleanResultsVector(max_ratio_algo, files);
    cleanResultsVector(min_ratio_algo, files);

}

void cleanResultsVector(std::vector <std::pair <std::pair <std::string , double > , int >> & max_ratio_algo, std::vector<std::string> files) {
    string algo, cleanAlgo;
    int points,cleanPoints;
    double  score_max, score_min, clean_score_max, clean_score_min;
    pair <std::string, double> pMax, pMin, cleanMax, cleanMin;
    int counter;
    vector <pair <pair <string , double > , int >>  new_max_ratio_algo;
    vector <pair <pair <string , double > , int >>  new_min_ratio_algo;
    for (int i = 0; i < max_ratio_algo.size(); i++) {
        pMax = max_ratio_algo[i].first; 
        points = max_ratio_algo[i].second; 
        algo = pMax.first;
        score_max = pMax.second;
        counter = 0; 
        for (int j = 0; j < new_max_ratio_algo.size(); j++) {
            cleanMax = new_max_ratio_algo[j].first; 
            cleanPoints = new_max_ratio_algo[j].second; 
            cleanAlgo = cleanMax.first;
            clean_score_max = cleanMax.second;
            if ((cleanAlgo == algo) && (clean_score_max == score_max) && (cleanPoints == points)) {
                counter++; 
            }
        }
        if (counter == 0) {
            new_max_ratio_algo.push_back(make_pair(pMax,points));
        }
    }
    
    max_ratio_algo = new_max_ratio_algo; 
}

double findMaxScore(std::vector <std::pair <std::pair <std::string , double > , int >> & max_ratio_algo, int points, string algorithm) {
    double max = 0.0;
    for (int i = 0; i < max_ratio_algo.size(); i++) {
        pair <std::string , double> p;
        p = max_ratio_algo[i].first; 
        string algo = p.first;
        double score = p.second; 
        int filepoints = max_ratio_algo[i].second;
        if ((algo == algorithm) && (points == filepoints)) {
            if (score > max) {
                max = score; 
            }
        }
    }
    return max; 
}

double findMinScore(std::vector <std::pair <std::pair <std::string , double > , int >> & min_ratio_algo, int points, string algorithm) {
    double min = 100.0;
    for (int i = 0; i < min_ratio_algo.size(); i++) {
        pair <std::string , double> p;
        p = min_ratio_algo[i].first; 
        string algo = p.first;
        double score = p.second; 
        int filepoints = min_ratio_algo[i].second;
        if ((algo == algorithm) && (points == filepoints)) {
            if (score < min) {
                min = score; 
            }
        }
    }
    return min; 
}

void printResultsBoard(std::vector <std::pair <std::pair <std::string , double > , int >> & max_ratio_algo,
                    std::vector <std::pair <std::pair <std::string , double > , int >> & min_ratio_algo, std::vector<std::string> files,
                    std::vector <std::pair <std::pair <std::string, double>, int>> & MAX_bounds, std::vector <std::pair <std::pair <std::string, double>, int>> & MIN_bounds)
{   
    
    sort(max_ratio_algo.begin(), max_ratio_algo.end(), sortbysec);
    int points; 
    string algo;
    double score_max; 
    vector <string> algo_results; 
    vector<int> algo_points; 
    std::vector<string>::iterator it;
    std::vector<int>::iterator points_it;

    for (int i = 0; i < max_ratio_algo.size(); i++) {
        pair <std::string, double> p;
        p = max_ratio_algo[i].first; 
        points = max_ratio_algo[i].second;
        algo = p.first;
        score_max = p.second;
        it = find(algo_results.begin(), algo_results.end(), algo);
        if (it != algo_results.end()) { 
            continue;
        }
        else { 
            algo_results.push_back(algo); 
        }
    }

    for (int i = 0; i < files.size(); i++) {
        int fp = getSumPoints(files[i]); 
        points_it = find(algo_points.begin(), algo_points.end(), fp);
        if (points_it != algo_points.end()) { 
            continue;
        }
        else { 
            algo_points.push_back(fp); 
        }
    }
    sort(algo_points.begin(), algo_points.end());
    
    
    int filepoints;
    double scoreMax,scoreMin,boundMax,boundMin;
    for (int i = 0; i < algo_results.size(); i++) {
        cout << "\t\t";
        cout << algo_results[i];
        cout << endl; 
        cout << "Size \t"; 
        cout << "||"; 
        cout << "min score | max_score | min_bound | max_bound ||" << endl;   
        for (int j = 0; j < algo_points.size(); j++) {
            scoreMax = findMaxScore(max_ratio_algo, algo_points[j], algo_results[i]); 
            scoreMin = findMinScore(min_ratio_algo, algo_points[j], algo_results[i]); 
            boundMax = findBoundMax(MAX_bounds, algo_points[j], algo_results[i]); 
            boundMin = findBoundMin(MIN_bounds, algo_points[j], algo_results[i]); 
            cout << algo_points[j] << "\t||" << scoreMin << "   |  " << scoreMax << "  |  " << boundMin << " | " << boundMax << endl; 
        }
        cout << endl; 
    }
    
}

bool sortbysec(const std::pair<std::pair<std::string,double>, int> &a,
              const std::pair<std::pair<std::string,double>, int> &b)
{
    return (a.second < b.second);
}


double findBoundMax(std::vector <std::pair <std::pair <std::string, double>, int>> & MAX_bounds, int points, string algorithm) {
    double max = 1.0;
    for (int i = 0; i < MAX_bounds.size(); i++) {
        pair <std::string , double> p;
        p = MAX_bounds[i].first; 
        if ((p.first == algorithm) && (MAX_bounds[i].second == points)) {
            if (p.second < max) {
                max = p.second; 
            }
        }
    }
    return max; 
}

double findBoundMin(std::vector <std::pair <std::pair <std::string, double>, int>> & MIN_bounds, int points, string algorithm) {
    double min = 0.0;
    for (int i = 0; i < MIN_bounds.size(); i++) {
        pair <std::string , double> p;
        p = MIN_bounds[i].first; 
        if ((p.first == algorithm) && (MIN_bounds[i].second == points)) {
            if (p.second > min) {
                min = p.second; 
            }
        }
    }
    return min; 
}

void createBoundVectors(std::vector <std::pair <std::string , int>> & algo_points, std::vector <std::pair <double, double>> & algo_res, 
                        std::vector <std::pair <std::pair <std::string, double>, int>> & MAX_bounds, std::vector <std::pair <std::pair <std::string, double>, int>> & MIN_bounds)
{
    for (int i = 0; i < algo_points.size(); i++) {
        for (int j = 0; j < algo_points.size(); j++) {
            if ((algo_points[i].first == algo_points[j].first) && (algo_points[i].second == algo_points[j].second)) {
                if (algo_res[j].first < algo_res[i].first) {
                    pair <std::string,double> p;
                    p.first = algo_points[i].first;
                    p.second = algo_res[j].first; 
                    int points = algo_points[i].second; 
                    MAX_bounds.push_back(make_pair(p,points)); 
                }
                if (algo_res[j].first > algo_res[i].first) {
                    pair <std::string,double> p;
                    p.first = algo_points[i].first;
                    p.second = algo_res[j].first; 
                    int points = algo_points[i].second; 
                    MIN_bounds.push_back(make_pair(p,points)); 
                }
            }
        }
    }
    sort(MAX_bounds.begin(), MAX_bounds.end()); 
    sort(MIN_bounds.begin(), MIN_bounds.end()); 
}


//simulated annealing 

bool is_intersected_P(Segment_2 ed, Polygon_2 A){   //if segment intersect polygon edges
    int count=0;
    for(const Segment_2& e: A.edges()){     //check if this segment intersect with the other segments of the polygon,(count>1)   
        if(CGAL::intersection(ed,e)) count++;  
    }
    if(count>0){
        return true;
    }else{
        return false;
    }
}

bool is_intersected_S(Segment_2 a, Segment_2 b){ //check if two segments are intersected      
    if(CGAL::intersection(a, b)){
        return true;
    }else{
        return false;
    }
}

// //function for finding the position of the polygon or set of points we want to insert or erase a point
int position(std::vector<std::pair<Point_2,int>> vertex_iterators, Point_2 x) {
   for(int i = 0; i<vertex_iterators.size(); i++){
        if ( x == vertex_iterators[i].first){
            return vertex_iterators[i].second; 
        }
    }
    return -1;
}

// //function for change the vertex iterators after insert/erase a point for a polygon
void change_after(std::vector<std::pair<Point_2,int>>& vertex_iterators ,Polygon_2 polygon) {
    vertex_iterators.clear(); 
    int i = 0;
    for (const Point_2& v: polygon.vertices()) {
        vertex_iterators.push_back(make_pair(v,i));
        i++; 
    }
}

//function to help us sort a polygon by x 
bool howtoSort(Point_2 p1, Point_2 p2) { return p1.x() < p2.x(); }

//finds the max x of a point od a polygon and returns the point
Point_2 find_maxx(Polygon_2 p){
    int maxx=p[0].x();
    Point_2 pos;
    for(int i=1;i<p.size();i++){
        if(p[i].x()>maxx) {
            maxx =p[i].x();
            pos = p[i];
        }
    }
    return pos;
}

//finds the min x of a point od a polygon and returns the point
Point_2 find_minx(Polygon_2 p){
    int minx=p[0].x();
    Point_2 pos;
    for(int i=1;i<p.size();i++){
        if(p[i].x()<minx) {
            minx =p[i].x();
            pos = p[i];
        }
    }
    return pos;
}

//function to compute the new polygon by simulated annealing
Polygon_2 simulated_annealing(Polygon_2 P,Polygon_2 ch, int n, string m, string an,int L, double cut_off, bool & exceed_cutoff){
    double T=1,E;
    Polygon_2 pl, max_polygon, min_polygon;

    for(auto it = P.begin(); it!= P.end();++it){
        pl.push_back(*it);  //pl is helping polygon where we initialize it with the elements of the P and make all the changes in it
    }

    std::vector<pair<Point_2,int>> vertex_iterators;
    for(int i=0; i<P.size(); i++){  //iterators for the polygonic line
        vertex_iterators.push_back(make_pair(P[i],i));
    }

    if(m=="-min"){    
        E = n*(abs(P.area())/abs(ch.area()));
    }else if(m=="-max"){
        E = n*(1-abs(P.area())/abs(ch.area()));
    }

    for(auto it = pl.begin(); it!= pl.end();++it){
        max_polygon.push_back(*it);
    }
    for(auto it = pl.begin(); it!= pl.end();++it){
        min_polygon.push_back(*it); 
    }

    if(an=="local"){
        bool end=false;
        int count=0;
        static double total_cutoff = 0.0;
        while(T>=0 && end==false){ 
            clock_t func_start, func_end;
            func_start = clock();
            Tree tree;
            Points points,result;
            bool acc=false; 
            Point_2 p,q,r,s;
            int random_num=0;
            if(L>n){
                random_selection( pl.begin()+1, pl.end()-1, n,std::back_inserter(result));
            }else{
                random_selection( pl.begin()+1, pl.end()-1, L,std::back_inserter(result));
            }
            q = result[random_num];
            while(acc==false && count<L && random_num<n){   
                for(auto it = P.begin(); it!= P.end();++it){
                    tree.insert(*it);
                }
                    
                int posq = position(vertex_iterators,q);
                if(posq!=pl.size()&&posq!=0){  
                    for(int i=0;i<pl.size();i++){
                        if(i==posq-1){
                            p=pl[i];
                        }
                        if(i==posq+1){
                            r=pl[i];
                        }
                    }
                    int posr = position(vertex_iterators,r);
                    for(int i=0;i<pl.size();i++){
                        if(i==posr+1){
                            s=pl[i];
                        }
                    }
                    //finding the maxx and maxy of this points 
                    points.push_back(p);
                    points.push_back(q);
                    points.push_back(r);
                    points.push_back(s);
                    Segment_2 pr = Segment_2(p,r);
                    Segment_2 qs = Segment_2(q,s);

                if(!is_intersected_S(pr,qs)){
                    //range of the rectangle
                    int minx,miny,maxx,maxy;
                    minx=points[0].x();
                    maxx=points[0].x();
                    miny=points[0].y();
                    maxy=points[0].y();
                    for(int i=1; i<points.size(); i++){
                        if(points[i].x()>maxx) maxx=points[i].x();
                        if(points[i].x()<minx) minx=points[i].x();
                        if(points[i].y()>maxy) maxy=points[i].y();
                        if(points[i].y()<miny) miny=points[i].y(); //change
                    }
                    //making our box with th minx,miny and maxx,maxy
                    Fuzzy_iso_box approximate_range(Point_2(minx,miny),Point_2(maxx,maxy),0.1);
                    tree.search( std::back_inserter( result ), approximate_range);;
                    Points clear_result;
                    for(int i=0; i<result.size(); i++){
                        int c=0;
                        for (int j=i+1; j<result.size(); j++){
                            if(result[i]==result[j]) c++;
                        }
                        if(c==0) clear_result.push_back(result[i]);
                    }
                        //points in the box
                    bool sn=false; //segments not intersected
                    for(int i=0;i<clear_result.size();i++){ //change
                        for(const Segment_2& e: pl.edges()){
                            if(e.source()==clear_result[i]||e.target()==clear_result[i]){
                                if((is_intersected_S(e,Segment_2(p,r)))==false&&(is_intersected_S(e,Segment_2(q,s))==false)){
                                    sn = true;
                                }
                            }
                        }
                    func_end = clock();
                    double func_time_taken = double(func_end - func_start) / double(CLOCKS_PER_SEC);
                    func_time_taken *=15;
                    total_cutoff+=func_time_taken;
                    if (total_cutoff >= cut_off) {                          //if optimazition exceeds cutoff time, terminate execution 
                        exceed_cutoff == true; 
                        if (m=="-max"){
                            return max_polygon;
                        }else {
                            return min_polygon; 
                        }
                    }   
                    }

                    if(sn==true){  //this means->no intersection;
                        int posr=position(vertex_iterators,r);
                        int posq=position(vertex_iterators,q);
                        pl.insert(pl.begin()+posr,q);
                        pl.insert(pl.begin()+posq,r);
                        change_after(vertex_iterators,pl);
                        Polygon_2 chn;
                        CGAL::convex_hull_2( pl.begin(), pl.end(), std::back_inserter(chn));

                        if(m=="-max"&& abs(pl.area())>abs(P.area())){
                           if(abs(pl.area())>abs(max_polygon.area())&&abs(pl.area())/abs(chn.area())<1){
                                max_polygon.clear();
                                for(auto it = pl.begin(); it!= pl.end();++it){
                                    max_polygon.push_back(*it);  
                                }
                            }
                        }else if(m=="-min"&& abs(pl.area())<abs(P.area())){
                            if(abs(pl.area())<abs(min_polygon.area())&&abs(pl.area())/abs(chn.area())<1){
                                min_polygon.clear();
                                for(auto it = pl.begin(); it!= pl.end();++it){
                                    min_polygon.push_back(*it);  
                                }
                            }
                        }else{
                            pl.erase(pl.begin()+posr);
                            pl.erase(pl.begin()+posq);
                            change_after(vertex_iterators,pl);
                        }
                        random_num++;
                        q = result[random_num];
                    }

                }

            }
            count++;
            
        }
        T = T-1/L;
        Polygon_2 ch1;
        CGAL::convex_hull_2( pl.begin(), pl.end(), std::back_inserter(ch1) );
        double E1= n*(abs(pl.area())/abs(ch1.area()));
        double de=E-E1;
        if(de<0){
            end = true;
            break;
            }else{
                double R=(float)rand() / (float)RAND_MAX; //generate r from 0->1
                if(exp(-de/T)>R) {  //metropolis
                    end = true;
                    break;
                }
            }  
        }
        if(m=="-max")  {
            return max_polygon;
        }else{
            return min_polygon;
        }
    }else if (an=="global"){
        bool end=false;
        static double total_cutoff = 0.0;
        while(T>=0 && end ==false){
            clock_t func_start, func_end;
            func_start = clock();
            int count=0;
            bool acc=false;
            int random_num=0;
            Points result;
            Point_2 p,q,r,s,t;
            if(L>n){
                random_selection( pl.begin()+1, pl.end()-1, n,std::back_inserter(result));
            }else{
                random_selection( pl.begin()+1, pl.end()-1, L,std::back_inserter(result));
            }
            q = result[random_num];
            s = result[random_num++];
            while(acc==false&&count<L){
                
                for(const Segment_2& e: pl.edges()){ //if qs is a polygon segment we are starting all over again
                    if(Segment_2(q,s) == e){
                         acc = true;
                    }
                }
                int posq = position(vertex_iterators,q);
                int poss = position(vertex_iterators,s);
                //finding p,r and q
                if(posq==0|| posq==pl.size()||poss==pl.size()) acc = true; //if q is at the first or at the end of pl we can;t find the p,r and if s is in the end we can't find t 
                    if(acc==false){
                    for(int i=0;i<pl.size();i++){
                        if(i==posq-1){
                            p=pl[i];
                        }
                        if(i==posq+1){
                            r=pl[i];
                        }
                        if(i==poss+1){
                            t=pl[i];
                        }
                    }
                    Segment_2 rp = Segment_2(r,p);
                    Segment_2 sq = Segment_2(s,q);
                    Segment_2 tq = Segment_2(t,q);
                    //checking if pr intersects sq or tq or any other segment, if its not we erase q from its position and we inserted it after s
                    if(is_intersected_S(rp,sq)==false&&is_intersected_S(rp,tq)==false){
                        std::cout <<"intersection1"<<std::endl;
                        if(is_intersected_P(rp,pl)==false){
                            std::cout <<"intersection2"<<std::endl;
                            // pl.erase(pl.begin()+posq);
                            pl.insert(pl.begin()+posq,r);
                            pl.insert(pl.begin()+poss,q);
                            change_after(vertex_iterators,pl);
                            Polygon_2 chn;
                            CGAL::convex_hull_2( pl.begin(), pl.end(), std::back_inserter(chn));
                                
                            if(m=="-max"&& abs(pl.area())>abs(P.area())){
                                if(abs(pl.area())>abs(max_polygon.area())&&abs(pl.area())/abs(chn.area())<1){
                                    max_polygon.clear();
                                    for(auto it = pl.begin(); it!= pl.end();++it){
                                        max_polygon.push_back(*it);  
                                    }
                                }
                            }else if(m=="-min"&& abs(pl.area())<abs(P.area())){
                                if(abs(pl.area())<abs(min_polygon.area())&&abs(pl.area())/abs(chn.area())<1){
                                    min_polygon.clear();
                                    for(auto it = pl.begin(); it!= pl.end();++it){
                                        min_polygon.push_back(*it);  
                                    }
                                }
                            }else{
                                pl.erase(pl.begin()+poss);
                                pl.erase(pl.begin()+posq);
                                change_after(vertex_iterators,pl);
                            }
                        }
                    }
                    random_num++;
                    q = result[random_num];
                    s = result[random_num++];
                    count++;
            }
            func_end = clock();
            double func_time_taken = double(func_end - func_start) / double(CLOCKS_PER_SEC);
            func_time_taken *=1000;
            total_cutoff+=func_time_taken;
            if (total_cutoff >= cut_off) {                          //if optimazition exceeds cutoff time, terminate execution 
                exceed_cutoff == true; 
                if (m=="-max"){
                    return max_polygon;
                }else {
                    return min_polygon; 
                }
                  
            }
        }
        T = T -1/L;
        Polygon_2 ch1;
        CGAL::convex_hull_2( pl.begin(), pl.end(), std::back_inserter(ch1) );
        double E1= n*(pl.area()/ch1.area());
        double de=E-E1;
        if(de<0){
            end = true;
            break;
        }else{
            double R=(float)rand() / (float)RAND_MAX; //generate r from 0->1
            if(exp(-de/T)>R) {  //metropolis
                end = true;
                break;
            }
        }  
    }
    if(m=="-max")  {
        return max_polygon;
    }else{
        return min_polygon;
    }
   }else if (an=="subdivisor"){
    //sorting our polygon
    int count=0;
    bool acc=false;
    static double total_cutoff = 0.0;
    while(acc==false && count <L){
        clock_t func_start, func_end;
        func_start = clock();
        std::sort(pl.begin(), pl.end(), howtoSort);
        srand((unsigned)time(0));
        int M=rand()%100+10;
        int k=(n-1)/M-1;

        //making a vector of polygons and making it as an edge
        std::vector<Polygon_2> S;
        int a=0,i=0;
        
        while(i<=k && a<=pl.size()){
            Polygon_2 pol;
            for(int j=0;j<M;j++){
                if(i==0){
                    pol[j] = pl[a];
                }else{
                    if(j==0){
                        pol[j] = pl[a-1]; //the last point of the previous set
                    }else{
                        pol[j] =pl[a];
                    }
                }
                a++;
                if(j==M-1 && a<pl.size()){  //vazoume ta upoloipa simeia sto teleutaio set
                    for(int l=a+1;l<pl.size(); l++){
                        pol[++j] = pl[l];

                    }
                }
            }
            for(int j=0; j<pol.size(); i++){
                S[i].push_back(pol[j]);
            }
            i++;
        }

         //finding the left and right segments of every set
        std::vector<Segment_2> rightsegs;
        std::vector<Segment_2> leftsegs;

        for(int i=0;i<k;i++){
            if(i!=0){
                Point_2 max=find_maxx(S[i]);
               for(const Segment_2& ed : S[i].edges()){
                    if(ed.target()==max) leftsegs.push_back(ed);
               }
            }
            if(i==k-1){
                Point_2 min=find_minx(S[i]);
                for(const Segment_2& ed : S[i].edges()){
                    if(ed.target()==min ) rightsegs.push_back(ed);
                }
            }
        }
        //for every set finding the right and the left segment of the next set and remove the r point
        Segment_2 pq,qr;
        for(int i=0;i<=k-1;i++){ 
            pq=leftsegs[i];
            qr=rightsegs[i+1];
            int posq=position(vertex_iterators,pq.target());
            int posp = position(vertex_iterators,pq.source());
            int posr = position(vertex_iterators,qr.target()); 
            pl.erase(pl.begin()+posr);
            // change_after(vertex_iterators,pl);
            if(pl.is_simple()){
                std::cout <<"im simple"<<std::endl;
                change_after(vertex_iterators,pl);
                acc=true;
            }else{
                pl.insert(pl.begin()+posr,qr.target());
                change_after(vertex_iterators,pl);
            }
        }
        count++;
        func_end = clock();
        double func_time_taken = double(func_end - func_start) / double(CLOCKS_PER_SEC);
        func_time_taken *=1000;
        total_cutoff+=func_time_taken;
        cout << "cutoff vs total cutoff " << cut_off << " , " << total_cutoff << endl; 
        if (total_cutoff >= cut_off) {                          //if optimazition exceeds cutoff time, terminate execution 
            return pl; 
        }
    }
    return pl;
    }else{
        std::cout<<"you should choose these three options: local,global,subdivisor"<<std::endl;
        return pl;
    } 
}