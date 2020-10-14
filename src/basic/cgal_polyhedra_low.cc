/*
Copyright (C) 2014 Adapted from CGAL example (Author: Pierre Alliez) by Vladimir Fuka
Copyright (C) 2020 Heiko Appel

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.
*/


#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;
typedef CGAL::Bbox_3 Bbox_3;

typedef struct {double x,y,z;} d3;
typedef struct {Polyhedron *poly; Tree *tree;} Polytree;

using std::cout;
using std::endl;

extern "C" int debuglevel;

extern "C" {

  void polyhedron_from_file (Polyhedron **poly, const char *fname, int verbose, int * const ierr){
    Polyhedron *polyhedron = new Polyhedron;

    std::ifstream in(fname);

    if (verbose) {cout << " Reading file " << fname << " " << endl;}

    try {
      in >> *polyhedron;
    }
    catch(...) {
      *ierr = 2;
      return;
    }

    if (verbose) {
      cout << " facets: " << polyhedron->size_of_facets() << endl;
      cout << " halfedges: " << polyhedron->size_of_halfedges() << endl;
      cout << " vertices: " << polyhedron->size_of_vertices() << endl;
    }

    if (polyhedron->size_of_facets()==0 ||
        polyhedron->size_of_halfedges()==0 ||
        polyhedron->size_of_vertices()==0){
          *ierr = 1;
          return;
        };

    *poly = polyhedron;
    *ierr = 0;
  }

  bool polyhedron_point_inside(Polyhedron *polyhedron, Point *query) {
    // Construct AABB tree with a KdTree
    Tree tree(faces(*polyhedron).first, faces(*polyhedron).second, *polyhedron);
    tree.accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    Point_inside inside_tester(tree);

    // Determine the side and return true if outside!
    return inside_tester(*query) == CGAL::ON_BOUNDED_SIDE || inside_tester(*query) == CGAL::ON_BOUNDARY;
  }

  void polyhedron_finalize(Polytree **pptree){
    delete (*pptree)->tree; (*pptree)->tree = NULL;
    delete (*pptree)->poly; (*pptree)->poly = NULL;
    delete *pptree; *pptree = NULL;
  }

}
