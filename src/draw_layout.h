/**
 * Declare layout drawing header
 */


#ifndef DRAW_LAYOUT_H
#define DRAW_LAYOUT_H




#include <stdlib.h>  // for exit()
#include <math.h>

#include <IL/ilut.h> // for output png


#include "config.h"
#include <vector>
#include <iostream>
#include <GLFW/glfw3.h>


/* Operation */
#define MAX_WIDTH_INIT_VAL  -10000000
#define BEZIER_INCREMENT    0.1

/* Property */
#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800


void draw_layout(std::vector< std::vector<VtxType> >& edges,
    std::vector< std::vector<CoordType> >& coord,
    std::vector<int>& clusters,
    std::vector< std::vector<CoordType> >& centers,
    std::vector<WgtType>& radius);


#endif