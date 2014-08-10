#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>

#include <iostream>
#include <iterator>
#include <fstream>
#include <list>
#include <cassert>
#include <map>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;

typedef CGAL::Fixed_alpha_shape_vertex_base_3<Gt>           Vb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<Gt>             Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>         Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>              Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>                Alpha_shape_3;
typedef CGAL::Fixed_alpha_shape_3<Triangulation_3>          Fixed_alpha_shape_3;

typedef Gt::Point_3                                         Point;
typedef Fixed_alpha_shape_3::Facet                          Facet;
typedef Fixed_alpha_shape_3::Cell_handle                    Cell_handle;

using namespace std;

int main()
{

    std::list<Point> lp;
    map<Point, int> pMap;

    // Read the input from stdin
    int nPts;
    double alpha;
    cin >> alpha;
    cerr << "Using alpha = " << alpha << "\n";
    cin >> nPts;
    cerr << "Reading " << nPts << " from stdin\n";
    Point p;
    int i = 0;
    for( ; nPts > 0; nPts--) {
        cin >> p; 
        lp.push_back(p);
        pMap[p] = i;
        i++;    
    }

    // Compute the alpha shape
    Fixed_alpha_shape_3 as(lp.begin(), lp.end(), alpha);
    //cerr << "Alpha shape computed in REGULARIZED mode by default\n";

    
    /*std::list<Facet> facetsReg;
    std::list<Facet> facetsExt;
    std::list<Facet> facetsInt;
    std::list<Facet> facetsSig;
    std::list<Facet> facetsAll;*/
    std::list<Facet> facetsOut;
    
    /*as.get_alpha_shape_facets(std::back_inserter(facetsReg), Fixed_alpha_shape_3::REGULAR);
    as.get_alpha_shape_facets(std::back_inserter(facetsExt), Fixed_alpha_shape_3::EXTERIOR);
    as.get_alpha_shape_facets(std::back_inserter(facetsInt), Fixed_alpha_shape_3::INTERIOR);
    as.get_alpha_shape_facets(std::back_inserter(facetsSig), Fixed_alpha_shape_3::SINGULAR);
    
    as.get_alpha_shape_facets(std::back_inserter(facetsAll), Fixed_alpha_shape_3::REGULAR);
    as.get_alpha_shape_facets(std::back_inserter(facetsAll), Fixed_alpha_shape_3::EXTERIOR);
    as.get_alpha_shape_facets(std::back_inserter(facetsAll), Fixed_alpha_shape_3::INTERIOR);
    as.get_alpha_shape_facets(std::back_inserter(facetsAll), Fixed_alpha_shape_3::SINGULAR);*/
    
    as.get_alpha_shape_facets(std::back_inserter(facetsOut), Fixed_alpha_shape_3::REGULAR);
    as.get_alpha_shape_facets(std::back_inserter(facetsOut), Fixed_alpha_shape_3::SINGULAR);
    
    /*cerr << "The alpha shape getter gives " << facetsReg.size() << " REGULAR facets\n";
    cerr << "The alpha shape getter gives " << facetsExt.size() << " EXTERIOR facets\n";
    cerr << "The alpha shape getter gives " << facetsInt.size() << " INTERIOR facets\n";
    cerr << "The alpha shape getter gives " << facetsSig.size() << " SINGULAR facets\n";
    cerr << "The alpha shape getter gives " << facetsAll.size() << " ALL facets\n";*/
    cerr << "The alpha shape getter gives " << facetsOut.size() << " OUT facets\n";
    
    // Print to stderr for debugging
    /*
    for(std::list<Facet>::const_iterator iter = facetsAll.begin(); iter != facetsAll.end(); iter++) {
        
        const Cell_handle& ch = iter->first;
        const int index = iter->second;

        const Point& a = ch->vertex((index+1)&3)->point();
        const Point& b = ch->vertex((index+2)&3)->point();
        const Point& c = ch->vertex((index+3)&3)->point();

        int cFace;
        cFace = as.classify(iter->first);
            
        cerr << pMap[a] << " ";
        cerr << pMap[b] << " ";
        cerr << pMap[c] << " ---- ";
        cerr << as.classify(iter->first) << "\n";
        
   }*/
    
   // Print to stdout for python program usage 
   for(std::list<Facet>::const_iterator iter = facetsOut.begin(); iter != facetsOut.end(); iter++) {
        
        const Cell_handle& ch = iter->first;
        const int index = iter->second;

        const Point& a = ch->vertex((index+1)&3)->point();
        const Point& b = ch->vertex((index+2)&3)->point();
        const Point& c = ch->vertex((index+3)&3)->point();
        
        
        cout << pMap[a] << " ";
        cout << pMap[b] << " ";
        cout << pMap[c] << "\n";
        
    }

    return 0;
}








