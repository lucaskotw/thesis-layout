#ifndef _INTRA_LAYOUT_H
#define _INTRA_LAYOUT_H


#include <vector>

#include "graph.h"
#include "config.h"
#include "eigenmat.h"

int intra_layout(
    Graph::Graph& g,
    int graphSize,
    DenseMat& distMat,
    int n_cls,
    std::vector<int>& clusters,
    double interpolation,
    std::vector< std::vector<CoordType> >& coord,
    std::vector< std::vector<CoordType> >& center_coord,
    std::vector< WgtType >& radii);


#endif