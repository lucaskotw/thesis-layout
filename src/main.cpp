#include <vector>
#include <iostream>
#include <fstream>
#include <regex>
#include <string>


#include "graph.h"
#include "load_graph.h"
#include "load_clusters.h"
#include "distance.h"
#include "intra_layout.h"
#include "inter_layout.h"
#include "draw_layout.h"

#define TEST_FLAG false

static
void output_layout(char* outfilePath,
    std::vector< std::vector<double> >& coord,
    std::vector< std::vector<double> >& center_coord,
    std::vector<double>& radii)
{
	using namespace std;

    int i; // index for loop

    string coord_outfile = string(outfilePath) + ".coord";
    string outdir = "out/";
    cout << coord_outfile << endl;

    // node coordinates
    fstream fo;
    fo.open(outdir+coord_outfile, fstream::out);
    for (i=0; i<coord.size(); ++i)
    {
        fo << coord.at(i).at(0) << " " << coord.at(i).at(1) << endl;
    }
    fo.close();

    // center coordinates
    string center_outfile = string(outfilePath) + ".center";
    fo.open(outdir+center_outfile, fstream::out);
    for (i=0; i<center_coord.size(); ++i)
    {
        fo << center_coord.at(i).at(0) << " " << center_coord.at(i).at(1) << endl;
    }
    fo.close();

    // radius
    string radii_outfile = string(outfilePath) + ".radii";
    fo.open(outdir+radii_outfile, fstream::out);
    for (i=0; i<radii.size(); ++i)
    {
        fo << radii.at(i) << endl;
    }
    fo.close();

}


int main(int argc, char ** argv)
{
	using namespace std;

    /* Load the Graph*/
    Graph g(0);
    // load_graph_from_gml(argv[1], g);
    load_graph_from_mm(argv[1], g);
    g.print_graph();

    // create edges
    vector< vector<int> > edges;
    vector<int> nbors;
    vector<int> edge(2, 0);
    for (int i=0; i<g.get_num_vtxs(); i++)
    {
        nbors = g.adj(i);
        for (int nb=0; nb<nbors.size(); ++nb)
        {
            if (i<nbors.at(nb))
            {
                edge.at(0) = i;
                edge.at(1) = nbors.at(nb);
                edges.push_back(edge);
            }
        }
    }
    cout << "edges size=" << edges.size() << endl;

    // create distance matrix
    DenseMat dist_mat(g.get_num_vtxs(), g.get_num_vtxs());
    distance_matrix(g, dist_mat);
    cout << "distance matrix" << endl;
    cout << dist_mat << endl;
    vector< vector<double> > coord(g.get_num_vtxs(), vector<double>(2));
    for (int c=0; c<2; c++)
    {
        for (int r=1; r<coord.size(); ++r)
        {
            coord.at(r).at(c) = rand()%100/50.0;
        }
    }


    // Clusters the Graph 
    vector<int> clusters(g.get_num_vtxs());
    int n_cls;
    load_clusters_from_group(argv[2], clusters, n_cls);

    cout << "cluster size = " << n_cls << endl;
    cout << "cluster vector size = " << clusters.size() << endl;
    cout << "cluster" << endl;
    for (vector<int>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    {
        cout << *it << " ";
    }
    cout << endl;


    /* Layout Algorithm */
    vector< vector<CoordType> > center_coord(n_cls, vector<CoordType>(2));

    vector< WgtType > radii(n_cls);

    cout << "center size=" << center_coord.size() << endl;
    cout << "radii size=" << radii.size() << endl;
    // double i_p = 1.0 / (double)(g.get_num_vtxs()+1);  // interpolation coeff
    const double i_p = 0.3;  // interpolation coeff

    intra_layout(g, g.get_num_vtxs(),
                 dist_mat, n_cls, clusters,
                 i_p, coord, center_coord, radii);
    inter_layout(argv[3], g, dist_mat, center_coord, radii, n_cls, clusters, coord); 


    /* Output the information and coordinates */
    cout << "coordinates" << endl;
    for (int i=0; i<coord.size(); ++i)
    {
        cout << coord.at(i).at(0) << " " << coord.at(i).at(1) << endl;
    }
    cout << "edges" << endl;
    for (int i=0; i<edges.size(); ++i)
    {
        cout << edges.at(i).at(0) << " " << edges.at(i).at(1) << endl;
    }
    cout << "center_coord" << endl;
    for (int i=0; i<center_coord.size(); ++i)
    {
        cout << center_coord.at(i).at(0) << " " << center_coord.at(i).at(1) << endl;
    }
    cout << "radii" << endl;
    for (int i=0; i<radii.size(); ++i)
    {
        cout << radii.at(i) << endl;
    }

    output_layout(argv[3], coord, center_coord, radii);
    draw_layout(edges, coord, clusters, center_coord, radii);

    



	return 0;
}