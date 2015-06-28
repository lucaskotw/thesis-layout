#include <vector>
#include <iostream>


#include "graph.h"
#include "load_graph.h"
#include "distance.h"

int main(int argc, char ** argv)
{
	using namespace std;

    /* Load the Graph*/
    Graph g(0);
    // load_graph_from_gml(argv[1], g);
    load_graph_from_mm(argv[1], g);
    g.print_graph();

    // create distance matrix
    DenseMat dist_mat(g.get_num_vtxs(), g.get_num_vtxs());
    distance_matrix(g, dist_mat);



	return 0;
}