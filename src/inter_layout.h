#ifndef _INTER_LAYOUT_H
#define _INTER_LAYOUT_H

#include "graph.h"
#include "config.h"
#include "eigenmat.h"

#include <vector>
#include <utility>


void force_direct_placement(Graph::Graph& g,
	DenseMat& dist_mat,
	std::vector< std::pair<int, int> >& edges,
	std::vector< WgtType >& radii,
	std::vector< std::vector<CoordType> >& coord,
	int iteration,
	double initStep);

int inter_layout(char* outfilePath, Graph::Graph& g,
	DenseMat& distMat,
    std::vector< std::vector<CoordType> >& center_coord,
    std::vector< WgtType >& radii,
    int cls,
    std::vector<int>& clusters,
    std::vector< std::vector<CoordType> >& coord);

#endif