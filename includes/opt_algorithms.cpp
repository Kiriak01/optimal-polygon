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

void opt_local_search(string filename, string algorithm, int L, string area_polygonization, string threshold,
                      string init_algo, int init_edge_selection, string init_vertex_sort, string annealing) { 
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

    // string algorithm = "local_search"; 

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
    
    

    double initialPolygonArea = abs(initialPolygon.area());  
    Polygon_2 initialConvexHull = calc_convex_hull(initialPolygon); 
    double initialCHarea = abs(initialConvexHull.area());
    double initialRatio = initialPolygonArea/initialCHarea; 
  
    std::vector<pair<Point_2,int>> vertex_iterators;                //initializing a vector of pairs that holds each vertex of polygon and the 
    reform(vertex_iterators,initialPolygon);                         //correct iterator to acquire it.
 

    // int L = stoi(argv[8]); 
    // int L = 5; 
    // string area_polygonization = argv[9]; 
    // string area_polygonization = "-max"; 
    bool max_area_polygonization,min_area_polygonization; 
    string m; 
    if (area_polygonization == "-max") {
        max_area_polygonization = true; 
        min_area_polygonization = false; 
        m = "max"; 
    }else if (area_polygonization == "-min") {
        max_area_polygonization = false; 
        min_area_polygonization = true; 
        m = "min"; 
    }
    
    int total_modifications = 0; 

    cout << endl <<  "Testing for filename: " << filename << endl; 

    if (algorithm == "local_search") {
    clock_t start, end;
    start = clock();
    float DA;
    // string fs(argv[11]); 
    string fs = "0.1"; 
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


    // for(const Point_2& vertex: initialPolygon.vertices()) {
    //     cout << vertex << endl; 
    // }
    // for(const Segment_2& edge: initialPolygon.edges()) {
    //     cout << edge << endl; 
    // }

    cout << "\nL: " << L << endl; 
    cout << "Optimized polygonization (2nd assign): " << area_polygonization << endl; 
    cout << "Algorithm: " << algorithm << endl; 
    cout << "area_initial: " << initialPolygonArea << endl; 
    cout << "area: " << abs(initialPolygon.area()) << endl; 
    cout << "ratio_initial: "<< initialRatio << endl; 

    Polygon_2 finalCH = calc_convex_hull(initialPolygon); 
    double finalCHarea = abs(finalCH.area()); 
    cout << "ratio: " << abs(initialPolygon.area())/finalCHarea << endl; 
    cout << "Construction time: " << fixed << round(time_taken) << setprecision(5);
    cout << " ms " << endl;    
    cout << "Initial polygon made with algorithm: " << init_algo << endl; 
    } else {

        clock_t start, end;
        start = clock();
        Polygon_2 ch;
        CGAL::convex_hull_2( initialPolygon.begin(), initialPolygon.end(), std::back_inserter(ch));
        int n=initialPolygon.size();
        // string ann = "subdivisor"; // 
        // string m = "min";
        Polygon_2 new_p = simulated_annealing(initialPolygon,ch,n,m,annealing,L); 
        end = clock();
        
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        time_taken *=1000;

        Polygon_2 chn;
        CGAL::convex_hull_2( initialPolygon.begin(), initialPolygon.end(), std::back_inserter(chn));
        
        cout << "Optimized polygonization (2nd assign): " << area_polygonization << endl;
        std::cout << "Algorithm: "<< algorithm << std::endl;
        std::cout << "area initial: " << initialPolygon.area() << std::endl;
        std::cout << "area: " << new_p.area() << std::endl;
        std::cout << "ratio_initial: " << initialPolygon.area()/ch.area() << std::endl;
        std::cout << "ratio: " << new_p.area()/chn.area() << std::endl;
        cout << "Construction time: " << fixed << round(time_taken) << setprecision(5);
        cout << " ms " << endl; 
        cout << "Initial polygon made with algorithm: " << init_algo << endl; 

    }
}; 




//simulated annealing 

bool is_intersected_P(Segment_2 ed, Polygon_2 A){   //if segment intersect polygon edges
    int count = 0;
    for(const Segment_2& e: A.edges()){     //check if this segment intersect with the other segments of the polygon,(count>1)   
        if(CGAL::do_intersect(ed,e)) count++;  
    }
    if(count >= 1){
        return true;
    }else{
        return false;
    }
}

bool is_intersected_S(Segment_2 a, Segment_2 b){ //check if two segments are intersected      
    if(CGAL::do_intersect(a, b)){
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
Polygon_2 simulated_annealing(Polygon_2 P,Polygon_2 ch, int n, string m, string an,int L){
    double T=1,E;
    int maxit=n;
    Polygon_2 pl;
    if(m=="min"){    
        E = n*(P.area()/ch.area());
    }else{
        E = n*(1-(P.area()/ch.area()));
    }

   for(auto it = P.begin(); it!= P.end();++it){
        pl.push_back(*it);  //pl is helping polygon where we initialize it with the elements of the P and make all the changes in it
    }

   std::vector<pair<Point_2,int>> vertex_iterators;
    for(int i=0; i<P.size(); i++){  //iterators for the polygonic line
      vertex_iterators.push_back(make_pair(P[i],i));
    }

    if(an=="local"){
        int count=0;
        while(T>=0 ){
            Tree tree;
            Points points,result;
            bool acc=false; //when this becomes true we don;t make any other iterations
            while(acc==false){
                for(auto it = P.begin(); it!= P.end();++it){
                    tree.insert(*it);
                }
                Point_2 p,q,r,s;
                random_selection( pl.begin(), pl.end(), 1,std::back_inserter(result));
                q = result[0];
                int posq = position(vertex_iterators,q);
            
                //finding p,r and s
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

               if(!is_intersected_S(Segment_2(p,r),Segment_2(r,s))){
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
                        if(points[i].y()<miny) maxx=points[i].y();
                    }

                    //making our box with th minx,miny and maxx,maxy
                    Fuzzy_iso_box range(Point_2(minx,miny),Point_2(maxx,maxy));
                    tree.search( std::back_inserter( result ), range);
                    //points in the box
                    int sn=0; //segments not intersected
                    for(int i=0;i<=result.size();i++){
                        for(const Segment_2& e: pl.edges()){
                            if(e.source()==result[i]||e.target()==result[i]){
                                if(!is_intersected_S(e,Segment_2(p,r))&&!is_intersected_S(e,Segment_2(q,s))){
                                    sn++;
                                }
                            }
                        }
                    }
                    if(sn==result.size()){  //this means->no intersections
                        int posr=position(vertex_iterators,r);
                        int posq=position(vertex_iterators,q);
                        pl.insert(pl.begin()+posr,q);
                        pl.insert(pl.begin()+posq,r);
                        change_after(vertex_iterators,pl);
                        acc=true;
                    }else{
                        L++;
                    }
                }
            }
        }
        T = T -1/L;
        Polygon_2 ch1;
        CGAL::convex_hull_2( pl.begin(), pl.end(), std::back_inserter(ch1) );
        double E1= n*(pl.area()/ch1.area());
        double de=E-E1;
        if(de<0){
            return pl;
        }else{
            double R=(float)rand() / (float)RAND_MAX; //generate r from 0->1
            if(exp(-de/T)>R) {  //metropolis
                return pl;
            }else{
               Polygon_2 p_end = simulated_annealing(pl,ch1,n,m,an,L);
            }
       }
     }
    }else if (an=="global"){
        int L=1;
        while(T>=0){
            bool acc=false;
            while(acc==false){
                Points points;
                Point_2 p,q,r,s,t;
                random_selection( pl.begin(), pl.end(), 2,std::back_inserter(points));
                q = points[0];
                s = points[1];
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
                    //checking if pr intersects sq or tq or any other segment, if its not we erase q from its position and we inserted it after s
                    if(!is_intersected_S(Segment_2(r,p),Segment_2(s,q))&&!is_intersected_S(Segment_2(r,p),Segment_2(t,q))){
                        if(!is_intersected_P(Segment_2(r,p),pl)){
                            pl.erase(pl.begin()+posq);
                            pl.insert(pl.begin()+posq,r);
                            pl.insert(pl.begin()+poss+1,q);
                            acc == true;
                        }
                    }else{
                        L++;
                    }
            }
        }
        T = T -1/L;
        Polygon_2 ch1;
        CGAL::convex_hull_2( pl.begin(), pl.end(), std::back_inserter(ch1) );
        double E1= n*(pl.area()/ch1.area());
        double de=E-E1;
        if(de<0){
            return pl;
        }else{
            double R=(float)rand() / (float)RAND_MAX; //generate r from 0->1
            if(exp(-de/T)>R) {  //metropolis
                return pl;
            }else{
               Polygon_2 p_end = simulated_annealing(pl,ch1,n,m,an,0);
            }
        }
    }
   }else if (an=="subdivisor"){
    //sorting our polygon
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
                change_after(vertex_iterators,pl);
            }
            return pl;
    }else{
        std::cout<<"you should choose these three options: local,global,subdivisor"<<std::endl;
        return pl;
    } 
}
