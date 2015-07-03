#ifndef _INTER_LAYOUT_H
#define _INTER_LAYOUT_H

#include "graph.h"
#include "config.h"
#include "eigenmat.h"

#include <vector>

int inter_layout(Graph::Graph& g,
	DenseMat& distMat,
    std::vector< std::vector<CoordType> >& center_coord,
    std::vector< WgtType >& radii,
    int cls,
    std::vector<int>& clusters,
    std::vector< std::vector<int> >& cluster_nodes_list,
    std::vector< std::vector<CoordType> >& coord);

#endif