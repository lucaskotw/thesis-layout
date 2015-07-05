#include "draw_layout.h"
#include <math.h>


#define PI 3.14159


/*******************
 * Basic Operation *
 *******************/
CoordType absolute(CoordType coordVal)
{
    if (coordVal < 0) return -coordVal;
    else return coordVal;
}


static
double get_drawing_width(std::vector< std::vector<CoordType> >& centers,\
    std::vector<WgtType>& radius)
{

    using namespace std;
    double max_width;
    vector<WgtType> bound_candi;
    double theta;
    double bound;
    for (int i=0; i<centers.size(); ++i)
    {
        theta = atan2( centers.at(i).at(1), centers.at(i).at(0) );
        bound = centers.at(i).at(0)+ radius.at(i) * cos(theta);
        bound_candi.push_back(bound);
        bound = centers.at(i).at(1)+radius.at(i) * sin(theta);
        bound_candi.push_back(bound);
    }
    max_width = max(\
                  *max_element(bound_candi.begin(), bound_candi.end()),\
                  fabs(*min_element(bound_candi.begin(), bound_candi.end()))\
                );
    cout << *min_element(bound_candi.begin(), bound_candi.end());
    cout << "max width=" << max_width << endl;
    max_width += 1; // offset

    

    return max_width;
}




/******************
 * Basic Property *
 ******************/
void set_points_attributes()
{
  glEnable(GL_POINT_SMOOTH);
  glPointSize(10);
}





/******************
 * Drawing Method *
 ******************/
void draw_vertices(double bound, std::vector< std::vector<CoordType> >& coord, \
    std::vector<int>& clusters)
{
  double coord_max = bound;


  glBegin(GL_POINTS);

    for (int i=0; i<coord.size(); ++i) {
        switch(clusters.at(i))
        {
            case 0:
                glColor3f(1.0f, 0.0f, 0.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            case 1:
                glColor3f(0.0f, 1.0f, 0.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            case 2:
                glColor3f(0.0f, 0.0f, 1.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            case 3:
                glColor3f(0.0f, 0.0f, 0.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            default:
                glColor3f(0.5f, 0.5f, 0.5f);  
                break;
        }

        glVertex2f(coord[i][0]/coord_max, coord[i][1]/coord_max);
        // std::cout << coord[i][0]/coord_max << ", " << coord[i][1]/coord_max << std::endl;
    }

  glEnd();
}



void draw_line(CoordType p1_x, CoordType p1_y,\
               CoordType p2_x, CoordType p2_y)
{
    glBegin(GL_LINES);
      glColor3f(0.5f, 0.5f, 0.5f);
      glVertex2f(p1_x, p1_y);
      glVertex2f(p2_x, p2_y);
    glEnd();
    glFlush();
}


/*
 * draw the edges
 */
 static
void draw_edges(double bound, std::vector< std::vector<VtxType> >& edges,
    std::vector< std::vector<CoordType> >& coord)
{
  double coord_max = bound;

  double pt1;
  double pt2;
  int cnt = -1; // count current dealing port and boundary id
  for (int i=0; i<edges.size(); ++i)
  {
    pt1 = edges.at(i).at(0);
    pt2 = edges.at(i).at(1);

    draw_line(coord.at(pt1).at(0)/coord_max, coord.at(pt1).at(1)/coord_max,
              coord.at(pt2).at(0)/coord_max, coord.at(pt2).at(1)/coord_max);
  }

  


}


void draw_radius(double bound, std::vector< std::vector<CoordType> >& coord,
                 std::vector< std::vector<CoordType> >& centers,
                 std::vector<WgtType>& radius)
{
  double coord_max = bound;
  double r;
  for (int c=0; c<centers.size(); ++c)
  {
    // draw radius
    r = radius.at(c);
    glBegin(GL_LINE_LOOP);
      switch(c)
        {
            case 0:
                glColor3f(1.0f, 0.0f, 0.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            case 1:
                glColor3f(0.0f, 1.0f, 0.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            case 2:
                glColor3f(0.0f, 0.0f, 1.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            case 3:
                glColor3f(0.0f, 0.0f, 0.0f);
                // glColor3f(1.0f, 1.0f, 1.0f);
                break;
            default:
                glColor3f(0.5f, 0.5f, 0.5f);  
                break;
        }
      for (int i=0; i<360; i++)
      {
        float degInRad = i*M_PI/180;
        glVertex2f((centers.at(c).at(0)+cos(degInRad)*r)/coord_max,\
                   (centers.at(c).at(1)+sin(degInRad)*r)/coord_max);
      }

    glEnd();

  }
}


void draw_layout(std::vector< std::vector<VtxType> >& edges,
    std::vector< std::vector<CoordType> >& coord,
    std::vector<int>& clusters,
    std::vector< std::vector<CoordType> >& centers,
    std::vector<WgtType>& radius)
{

  using namespace std;
  GLFWwindow* window;

  /* Initialize the library */
  if (!glfwInit())
      exit(1);
      // return -1;

  /* Initalize DevIL */
  // ilutRenderer(ILUT_OPENGL);

  /* Create a windowed mode window and its OpenGL context */
  window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "2D Stress Layout", NULL, NULL);
  if (!window)
  {
      glfwTerminate();
      exit(1);
      // return -1;
  }

  /* Make the window's context current */
  glfwMakeContextCurrent(window);


  /* Loop until the user closes the window */

  double screen_bounding_width = get_drawing_width(centers, radius);
  cout << "screen_bounding_width = " << screen_bounding_width << endl;

  while (!glfwWindowShouldClose(window))
  {
    /* Render here */

    // glfwGetFramebufferSize(window, &width, &height);
    // glViewport(0, 0, width, height);
    set_points_attributes();

    // set frame background color
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear(GL_COLOR_BUFFER_BIT);

#if 0
  draw_edges(screen_bounding_width, edges, centers);
    // draw_vertices(screen_bounding_width, coord, clusters);

    // for seeing the radius of clusters
    draw_radius(screen_bounding_width, coord, centers, radius);
#else
    draw_edges(screen_bounding_width, edges, coord);
    draw_vertices(screen_bounding_width, coord, clusters);

    // for seeing the radius of clusters
    draw_radius(screen_bounding_width, coord, centers, radius);
#endif

    /* Swap front and back buffers */
    glfwSwapBuffers(window);

    /* Poll for and process events */
    glfwPollEvents();
  }

  glfwTerminate();
}