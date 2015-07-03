#include "intra_layout.h"
#include "sm.h"



static
void create_cluster_dist_mat(
    DenseMat& distMat,
    int cls,
    std::vector< std::vector<int> >& cluster_nodes_list,
    std::vector<WgtType>& nodes_radii,
    DenseMat& clusterDistMat)
{
    // Steps
    // 1. resize cluster distance matrix
    // 2. match the distance from distMat to clusterDistMat
    // 3. add central node influence

    // Step 1
    clusterDistMat.resize(cluster_nodes_list.at(cls).size()+1, cluster_nodes_list.at(cls).size()+1);
    clusterDistMat.fill(0);

    // Step 2
    int rowIdx;
    int colIdx;
    for (int c=0; c<cluster_nodes_list.at(cls).size(); ++c)
    {
        for (int r=0; r<cluster_nodes_list.at(cls).size(); ++r)
        {
            rowIdx = cluster_nodes_list.at(cls).at(r);
            colIdx = cluster_nodes_list.at(cls).at(c);
            clusterDistMat(r+1, c+1) = distMat(rowIdx, colIdx);
        }
    }

    // Step 3
    for (int i=0; i<nodes_radii.size(); ++i)
    {
        clusterDistMat(0, i+1) = nodes_radii.at(i);
        clusterDistMat(i+1, 0) = nodes_radii.at(i);

    }
}


static
void match_cluster_coord(int cls, std::vector< std::vector<int> >& cluster_nodes_list,
                        std::vector< std::vector<CoordType> >& intra_coord,
                        std::vector< std::vector<CoordType> >& coord)
{
    // Assign intro_coord to coord based on cluster_vtxs
    int vtx_id;
    for (int i=0; i<cluster_nodes_list.at(cls).size(); ++i)
    {
        vtx_id = cluster_nodes_list.at(cls).at(i);
        coord.at(vtx_id).at(0) = intra_coord.at(i+1).at(0);
        coord.at(vtx_id).at(1) = intra_coord.at(i+1).at(1);
    }
}


static
void calculate_nodes_radii(Graph::Graph& g,
    DenseMat& distMat,
    int cls,
    std::vector<int>& clusters,
    std::vector< std::vector<int> >& cluster_nodes_list,
    std::vector< std::vector<CoordType> >& intra_coord,
    std::vector<WgtType>& nodes_radii)
{
    // Steps
    // 1. Calculate cluster radius
    // 2. get the inter/intra cluster degree of each vtxs
    // 3. calculate radial constriants
    using namespace std;

    // Step 1
    // get the corresponding cluster distance
    // minus central node
    DenseMat clsDistMat(intra_coord.size()-1, intra_coord.size()-1);
    VtxType rr;
    VtxType cc;
    for (int c=0; c<clsDistMat.cols(); ++c)
    {
        for (int r=0; r<clsDistMat.rows(); ++r)
        {
            rr = cluster_nodes_list.at(cls).at(r);
            cc = cluster_nodes_list.at(cls).at(c);
            clsDistMat(c, r) = distMat(cc, rr);
        }
    }
    cout << "cls distance matrix" << endl;
    cout << clsDistMat << endl;
    // find the maximum pair distance
    double cls_radius = clsDistMat.maxCoeff()/2;

    // Step 2
    vector<VtxType> intra_deg(cluster_nodes_list.at(cls).size(), 0);
    vector<VtxType> inter_deg(cluster_nodes_list.at(cls).size(), 0);
    vector<VtxType> nbors;
    for (int i=0; i<cluster_nodes_list.at(cls).size(); ++i)
    {
        nbors = g.adj( cluster_nodes_list.at(cls).at(i) );
        for (int n=0; n<nbors.size(); ++n)
        {
            if (clusters.at( cluster_nodes_list.at(cls).at(i) ) == clusters.at( nbors.at(n) ))
            {
                intra_deg.at(i) += 1;
            }
            else
            {
                inter_deg.at(i) += 1;
            }
        }
    }
    cout << "intra_deg=" << intra_deg.size() << endl;
    cout << "inter_deg=" << inter_deg.size() << endl;


    // Step 3
    const int offset = 1;
    double min_inter = *min_element(inter_deg.begin(), inter_deg.end());
    double max_inter = *max_element(inter_deg.begin(), inter_deg.end());
    cout << "radius" << endl;
    for (int i=0; i<nodes_radii.size(); ++i)
    {
        nodes_radii.at(i) = cls_radius*( (inter_deg.at(i)-min_inter+offset) / (max_inter-min_inter+offset) );
    }



}


static
WgtType calculate_norm(std::vector< CoordType >& pt1, std::vector< CoordType >& pt2)
{
	using namespace std;
    // process 2-norm of given input
    WgtType res;
    res = pow(pt1.at(0)- pt2.at(0), 2) + pow(pt1.at(1)- pt2.at(1), 2);

    return sqrt(res);
}


static
int get_radii_centers(
	int cls,
	std::vector< std::vector<CoordType> >& intra_coord,
    std::vector< WgtType >& radii,
    std::vector< std::vector<CoordType> >& center_coord)
{
    // Steps
    //   1. get center and get the radii
    //   2. assign to corresponds radii and center_coord
    std::vector<CoordType> center(2, 0);
    WgtType radius = 0; // radius default value
    WgtType candi_radius = 0; // candidate radius value


    for (int i=0; i<intra_coord.size(); ++i)
    {
    	// step 1
    	center.at(0) += intra_coord.at(i).at(0);
    	center.at(1) += intra_coord.at(i).at(1);

    	// step 2
        candi_radius = calculate_norm(intra_coord.at(i), center);
        if (candi_radius > radius) radius = candi_radius;
    }

    center.at(0) /= intra_coord.size()-1;
    center.at(1) /= intra_coord.size()-1;

    center_coord.at(cls) = center;
    radii.at(cls) = radius;


    return SUCCESS_GET_RADII_CENTERS;
}



int intra_layout(
    Graph::Graph& g,
    int graphSize,
    DenseMat& distMat,
    std::vector<int>& clusters,
    std::vector< std::vector<int> >& cluster_nodes_list,
    double interpolation,
    std::vector< std::vector<CoordType> >& coord,
    std::vector< std::vector<CoordType> >& center_coord,
    std::vector< WgtType >& radii)
{
    // Steps (with central node added)
    // Forloop from 0 to nCluster-1
    //     1. get the nodes inside cluster: O(V)
    //     2. layout inside cluster:
    //       - calculate vertices distance
    //       - create distance matrix
    //       - run stress majorization
    //     3. put the coordinates to the coord: O(N)
    using namespace std;
    vector< vector<CoordType> > intra_coord(0, vector<CoordType>(2));
    cout << "interpolation=" << interpolation << endl;

    vector<WgtType> nodes_radii;
    DenseMat cluster_dist_mat;
    
    for (int c=0; c<cluster_nodes_list.size(); ++c)
    {
        // loop step 1
        intra_coord.resize(cluster_nodes_list.at(c).size()+1, vector<CoordType>(2, 0));

        // random initialization
        intra_coord.at(0).at(0) = 0;
        intra_coord.at(0).at(1) = 0;
        for (int cc=0; cc<2; cc++)
        {
            for (int rr=1; rr<intra_coord.size(); ++rr)
            {
                intra_coord.at(rr).at(cc) = rand()%100/50.0;
            }
        }
        cout << "intra_coord assign complete" << endl;
        

        // loop step 2
        nodes_radii.resize(cluster_nodes_list.at(c).size());
        calculate_nodes_radii(g, distMat, c, clusters, cluster_nodes_list, intra_coord, nodes_radii);
        cout << "nodes radii size=" << nodes_radii.size() << endl;
        for (int i=0; i<nodes_radii.size(); ++i)
        {
            cout << nodes_radii.at(i) << " ";
        }
        cout << endl;
        std::cout << "finish calculate nodes radii" << std::endl;
        create_cluster_dist_mat(distMat, c, cluster_nodes_list, nodes_radii, cluster_dist_mat);
        cout << "cluster dist mat" << endl;
        cout << cluster_dist_mat << endl;
        stress_majorization_radial_refinement(cluster_dist_mat,
                                              interpolation,
                                              intra_coord);
        cout << "intra coord" << endl;
        for (vector< vector<CoordType> >::iterator it1=intra_coord.begin();\
            it1!=intra_coord.end();
            ++it1)
        {
            for (std::vector<CoordType>::iterator it2=(*it1).begin();\
            it2!=(*it1).end();
            ++it2)
            {
                std::cout << *it2 << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "cluster #" << c << " finish sm" << std::endl;

        // loop step 3
        match_cluster_coord(c, cluster_nodes_list, intra_coord, coord);

        // get radius and center
        get_radii_centers(c, intra_coord, radii, center_coord);
       
    }

    // print out center and radius
    cout << "clusters info" << endl;
    for (int i=0; i<center_coord.size(); ++i)
    {
    	cout << "x = " << center_coord.at(i).at(0) << ", y = " << center_coord.at(i).at(1)
    	<< ", r = " << radii.at(i) << endl;
    }
    
    

    return SUCCESS_INTRA_LAYOUT;
}