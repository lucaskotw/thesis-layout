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



static
void output_layout(char* outfilePath,
    std::vector< std::vector<double> >& coord,
    std::vector< std::vector<double> >& center_coord,
    std::vector<double>& radii)
{
	using namespace std;

	regex rgx(".*/(\\w+)/.*");
    smatch match;
    string data(outfilePath);
    regex_search(data, match, rgx);
    cout << match[1] << endl;
    string coord_outfile = string(match[1]) + ".coord";
    string outdir = "out/";
    cout << coord_outfile << endl;

    // node coordinates
    fstream fo;
    fo.open(outdir+coord_outfile, fstream::out);
    for (int i=0; i<coord.size(); ++i)
    {
        fo << coord.at(i).at(0) << " " << coord.at(i).at(1) << endl;
    }
    fo.close();

    // center coordinates
    string center_outfile = string(match[1]) + ".center";
    fo.open(outdir+center_outfile, fstream::out);
    for (int i=0; i<center_coord.size(); ++i)
    {
        fo << center_coord.at(i).at(0) << " " << center_coord.at(i).at(1) << endl;
    }
    fo.close();

    // radius
    string radii_outfile = string(match[1]) + ".radii";
    fo.open(outdir+radii_outfile, fstream::out);
    for (int i=0; i<radii.size(); ++i)
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

    // create distance matrix
    DenseMat dist_mat(g.get_num_vtxs(), g.get_num_vtxs());
    distance_matrix(g, dist_mat);
    vector< vector<double> > coord(g.get_num_vtxs(), vector<double>(2));
    for (int c=0; c<2; c++)
    {
        for (int r=1; r<coord.size(); ++r)
        {
            coord.at(r).at(c) = rand()%100/50.0;
        }
    }
#if 0
    // stress majorization
    

    stress_majorization(g.get_num_vtxs(), dist_mat, coord);


#else

    /* Clusters the Graph */
    std::vector<int> clusters(g.get_num_vtxs());
    int n_cls;
    load_clusters_from_group(argv[2], clusters, n_cls);

    cout << "cluster size = " << n_cls << endl;
    std::cout << "cluster" << std::endl;
    for (std::vector<int>::iterator it=clusters.begin();\
    it!=clusters.end();
    ++it)
    {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    // create cluster nodes list
    vector< vector<int> > cluster_nodes_list(n_cls, vector<int>(0));
	int cls;
    for (int i=0; i<clusters.size(); ++i)
    {
    	cls = clusters.at(i);
    	cout << cls << endl;
    	cluster_nodes_list.at(cls).push_back(i);
    }


	/* Layout Algorithm */
    vector< vector<CoordType> > center_coord(n_cls, vector<CoordType>(2));

    vector< WgtType > radii(n_cls);

    // double i_p = 1.0 / (double)(g.get_num_vtxs()+1);  // interpolation coeff
	const double i_p = 0.0;  // interpolation coeff

    intra_layout(g, g.get_num_vtxs(),
    			 dist_mat, clusters, cluster_nodes_list, i_p,
    			 coord, center_coord, radii);
    inter_layout(g, dist_mat, center_coord, radii, n_cls, clusters, cluster_nodes_list, coord);


#endif

    /* Output the information and coordinates */
    // node coordinates
    cout << "coord" << endl;
    for (int i=0; i<coord.size(); ++i)
    {
        cout << coord.at(i).at(0) << " " << coord.at(i).at(1) << endl;
    }

    cout << "center_coord" << endl;
    for (int i=0; i<center_coord.size(); ++i)
    {
        cout << center_coord.at(i).at(0) << " " << center_coord.at(i).at(1) << endl;
    }
	output_layout(argv[1], coord, center_coord, radii);



	return 0;
}