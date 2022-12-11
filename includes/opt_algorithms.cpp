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
        std::cout << "IN" <<std::endl;
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
