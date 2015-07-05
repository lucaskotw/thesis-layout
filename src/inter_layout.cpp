#include "inter_layout.h"
#include "bfs.h"
#include "distance.h"
#include "sm.h"
#include "log.h"
#include <math.h>
#include <algorithm>  // std::find, std::max_element
#include <utility>    // std::pair
#include <fstream>    // log out first few iteration coord, center, and radii
#include <string>     // string for output coord file
#include <boost/date_time/posix_time/posix_time.hpp>   // for current time check
#include <boost/date_time/posix_time/posix_time_io.hpp>


#define INFINITY_ENERGY 100000
#define IDEAL_BOUNDARY_FACTOR 0.1


/******************************************************************************
 *                      Global var in inter layout                            *
 ******************************************************************************/

static int inc_progress = 0;
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/******************************************************************************
 *                Standard library operator custom overload                   *
 ******************************************************************************/

// template <typename T,typename U>                                            
// std::pair<T,U> operator+(const std::pair<T,U> & u,const std::pair<T,U> & v)
// {
//     return {u.first+v.first,u.second+v.second};                                    
// }


// template <typename T,typename U>                                            
// std::pair<T,U> operator-(const std::pair<T,U> & u,const std::pair<T,U> & v)
// {
//     return {u.first-v.first, u.second-v.second};
// }


/******************************************************************************
 *                         Create Cluster Graph                               *
 ******************************************************************************/
static
void create_cluster_graph(Graph::Graph& g,
    std::vector<int> clusters,
    int nCluster,
    std::vector<WgtType>& radii,
    Graph::Graph& cg)
{
    // Steps
    //   1. get the linkage between cluster
    //   2. link the edge
    using namespace std;
    vector< vector<int> > cls_connection(nCluster, vector<int>(nCluster, 0));
    vector<int> cls_nodes;
    vector<int> nbors;


    // Step 1
    for (int c=0; c<nCluster; ++c)
    {
        cls_nodes.clear();
        cls_nodes.resize(0);

        for (int i=0; i<g.get_num_vtxs(); ++i)
            if (clusters.at(i) == c) cls_nodes.push_back(i);


        for (int i=0; i<cls_nodes.size(); ++i)
        {
            nbors = g.adj(cls_nodes.at(i));
            for (int n=0; n<nbors.size(); ++n)
            {
                if (clusters.at(nbors.at(n)) != c)
                    cls_connection.at(c).at(clusters.at(nbors.at(n))) += 1;
            }
        }
    }
    // Step 2
    for (int i=0; i<nCluster; ++i)
    {
        for (int j=i+1; j<nCluster; ++j)
            if (cls_connection.at(i).at(j) != 0) cg.add_edge(i, j, radii.at(i) + radii.at(j));
    }
}


static
void create_cluster_graph_distance(Graph::Graph& cg, DenseMat& cluster_dist_mat)
{

}


/******************************************************************************
 *                   Force-direct without node overlap                        *
 ******************************************************************************/
static
void update_step(double& step, double pre_energy, double curr_energy)
{
	if (curr_energy<pre_energy)
	{
		++inc_progress;
		if (inc_progress > 4)
		{
			inc_progress = 0;
			step /= 0.9;
		}
	}
	else
	{
		inc_progress = 0;
		step *= 0.9;
	}
}

static
double norm(std::vector<CoordType>& pt1, std::vector<CoordType>& pt2)
{
	using namespace std;
	double x_sqr = pow(pt1.at(0) - pt2.at(0), 2);
	double y_sqr = pow(pt1.at(1) - pt2.at(1), 2);
	return sqrt(x_sqr+y_sqr);
}

static
void calculate_energy(DenseMat& cluster_dist_mat,
	std::vector< std::vector<CoordType> >& center_coord,
	double& curr_energy)
{
	using namespace std;
	curr_energy = 0;

	vector<CoordType> u(2);
	vector<CoordType> v(2);
	for (int i=0; i<cluster_dist_mat.rows(); ++i)
	{
		u.at(0) = center_coord.at(i).at(0);
		u.at(1) = center_coord.at(i).at(1);
		for (int j=i+1; j<cluster_dist_mat.rows(); ++j)
		{
			v.at(0) = center_coord.at(j).at(0);
			v.at(0) = center_coord.at(j).at(1);
			curr_energy += pow(norm(u, v) - cluster_dist_mat(i, j) , 2);
		}
	}
}


static
double pos_diff_dist(std::pair<WgtType, WgtType>& pos_diff)
{
	using namespace std;
	return sqrt(pow(pos_diff.first, 2) + pow(pos_diff.second, 2));
}


static
double repulsive_force(std::pair<WgtType, WgtType>& pos_diff, double r_i, double r_j)
{
	using namespace std;
	// return pow(10+r_i+r_j, 2)/ pos_diff_dist(pos_diff);
	return pow((0.5)*(r_i+r_j), 2)/ pos_diff_dist(pos_diff);
}


static
double attractive_force(std::pair<WgtType, WgtType>& pos_diff, double r_i, double r_j)
{
	using namespace std;
	
	// return pow(pos_diff_dist(pos_diff), 2)/(10+r_i+r_j);
	return pow(pos_diff_dist(pos_diff), 2)/((0.5)*(r_i+r_j));
}


static
void calculate_rotate_angle(Graph::Graph& g, int cls, std::vector<int>& cluster_nodes,
	std::vector<int>& clusters,
	std::vector< std::vector<CoordType> >& coord,
	std::vector< std::vector<CoordType> >& center_coord,
	std::vector<WgtType>& radii,
	std::vector< WgtType >& rotate_angles)
{
	using namespace std;
	// torque
    int r_u; // current node id
    int r_v; // adjacent node id
    CoordType n_c_x; // center of neighbor's in different cluster
    CoordType n_c_y;
    CoordType r_u_x; // x components of raidus of u
    CoordType r_u_y;
    CoordType c_x;  // center of u's x coord
    CoordType c_y;
    vector<VtxType> nbors; // current node id
    double sin_coeff = 0.0;
    double cos_coeff = 0.0;
    pair<CoordType, CoordType> force;
    pair<CoordType, CoordType> arm;
    double force_val;
    double arm_val;
    double t_angle; // angle of torque

    // topo algo
    for (int j=0; j<cluster_nodes.size(); ++j)
    {
		r_u = cluster_nodes.at(j);
        r_u_x = coord.at(r_u).at(0);
        r_u_y = coord.at(r_u).at(1);
        c_x = center_coord.at(cls).at(0);
        c_y = center_coord.at(cls).at(1);
        arm = make_pair(r_u_x - c_x, r_u_y - c_y);
        arm_val = sqrt( pow(arm.first, 2) + pow(arm.second, 2));

		nbors = g.adj(r_u);
		for (int n=0; n<nbors.size(); ++n)
		{
			r_v = nbors.at(n);
			if ( clusters.at(r_u) != clusters.at(r_v) )
			{
				n_c_x = center_coord.at(clusters.at(r_v)).at(0);
            	n_c_y = center_coord.at(clusters.at(r_v)).at(1);

				// Rotate Step 1: calculate force and angle
				force = make_pair(n_c_x-r_u_x, n_c_y-r_u_y);
                force_val = sqrt( pow(force.first, 2) + pow(force.second, 2));
                force = make_pair(force.first/force_val, force.second/force_val);
                t_angle = (arm_val/radii.at(cls))*M_PI/2*sgn(arm.first*force.second-arm.second*force.first)*(arm.first*force.first+arm.second*force.second);
                rotate_angles.at(cls) += t_angle;
			}
		}

    }

	// for (int j=0; j<cluster_nodes.size(); ++j)
	// {
	// 	r_u = cluster_nodes.at(j);
 //        r_u_x = coord.at(r_u).at(0);
 //        r_u_y = coord.at(r_u).at(1);
 //        c_x = center_coord.at(cls).at(0);
 //        c_y = center_coord.at(cls).at(1);
 //        arm = make_pair(r_u_x - c_x, r_u_y - c_y);
 //        arm_val = sqrt( pow(arm.first, 2) + pow(arm.second, 2));

	// 	nbors = g.adj(r_u);
	// 	for (int n=0; n<nbors.size(); ++n)
	// 	{
	// 		r_v = nbors.at(n);
	// 		if ( clusters.at(r_u) != clusters.at(r_v) )
	// 		{
	// 			n_c_x = center_coord.at(clusters.at(r_v)).at(0);
 //            	n_c_y = center_coord.at(clusters.at(r_v)).at(1);

	// 			// Rotate Step 1: calculate force and angle
	// 			force = make_pair(n_c_x-r_u_x, n_c_y-r_u_y);
 //                force_val = sqrt( pow(force.first, 2) + pow(force.second, 2));
 //                t_angle = acos( (arm.first*force.first+arm.second*force.second) / (arm_val*force_val) );

 //                // Rotate Step 2
 //                // force_val = 1; // make force to be unit
 //                force_val = 1/force_val; // make force to be inverse to the distance
 //                sin_coeff += arm_val*force_val*cos(t_angle);
 //                cos_coeff += arm_val*force_val*sin(t_angle);
	// 		}
	// 	}
	// }

	// // Step 3
 //    cout << sin_coeff << " " << cos_coeff << endl;
 //    rotate_angles.at(cls) = atan(-cos_coeff/sin_coeff) * M_PI / 180;

}


static
void force_direct_with_torque(char* outfilePath, Graph::Graph& g, Graph::Graph& cg,
	std::vector< std::pair<int, int> > c_edges,
	DenseMat& cluster_dist_mat,
	std::vector<int>& clusters,
	std::vector< std::vector<CoordType> >& center_coord,
	std::vector< WgtType >& radii,
	std::vector< std::vector<CoordType> >& coord,
	int iteration,
	double initStep)
{

	// Iteration Steps
	// 1. calculate repulsive displacement
	// 2. calculate attractive displacement
	// 3. calculate rotate angle
	// 4. shift and rotate the cluster nodes and corresponding nodes
	using namespace std;
	using namespace boost::posix_time;

	// algorithm declaration
	double step = initStep;
	double curr_energy = INFINITY_ENERGY;
	double pre_energy;

	// for iteration dumping
	fstream fo;

	// [modified!] clusters_nodes
	vector<VtxType> cluster_nodes; // id: vtx id of clustered graph
                                  // val: vtx id in whole graph

	// force equilibrium declaration
	vector< pair<WgtType, WgtType> > node_disp(cg.get_num_vtxs());
	vector< WgtType > rotate_angles(cg.get_num_vtxs());
	pair<WgtType, WgtType> pos_diff;
	double rep_force;
	double att_force;
	double disp_x;
	double disp_y;
	double disp_val;
	int e_u;  // u in c-edge (u, v)
	int e_v;  // v in c-edge (u, v)
	double disp_u_x;
	double disp_u_y;
	double disp_v_x;
	double disp_v_y;
	

	CoordType c_x;  // center of u's x coord
    CoordType c_y;
    CoordType new_x; // relative to center nodes' x coordinates
    CoordType new_y; // relative to center nodes' x coordinates

    double angle;   // perform rotation for current node

    // dump before force process
    fo.open("out/"+string(outfilePath)+"_start.coord", fstream::out);
	for (int i=0; i<coord.size(); ++i)
	{
		fo << coord.at(i).at(0) << " " << coord.at(i).at(1) << endl;
	}
	fo.close();
	// center coord
	fo.open("out/"+string(outfilePath)+"_start.center", fstream::out);
	for (int i=0; i<center_coord.size(); ++i)
	{
		fo << center_coord.at(i).at(0) << " " << center_coord.at(i).at(1) << endl;
	}
	fo.close();

	// angle
	fo.open("out/"+string(outfilePath)+"_start.angle", fstream::out);
	for (int i=0; i<rotate_angles.size(); ++i)
	{
		fo << rotate_angles.at(i) << endl;
	}
	fo.close();

	for (int t=0; t<iteration; t++)
	{
		pre_energy = curr_energy;
		curr_energy = 0;
		// step 1: repulsive force placement
		for (int i=0; i<cg.get_num_vtxs(); ++i)
		{
			node_disp.at(i) = make_pair(0, 0);
			for (int j=0; j<cg.get_num_vtxs(); ++j)
			{
				if (i != j)
				{
					pos_diff = make_pair(center_coord.at(i).at(0)-center_coord.at(j).at(0),
										 center_coord.at(i).at(1)-center_coord.at(j).at(1));
					rep_force = repulsive_force(pos_diff, radii.at(i), radii.at(j));
					disp_x = node_disp.at(i).first + pos_diff.first/pos_diff_dist(pos_diff)*rep_force;
					disp_y = node_disp.at(i).second + pos_diff.second/pos_diff_dist(pos_diff)*rep_force;
					node_disp.at(i) = make_pair(disp_x+node_disp.at(i).first, disp_y+node_disp.at(i).second);
					if (t==0)
					{
						cout << "pos diff of cv_" << i << " = " << pos_diff.first << ", " << pos_diff.second << endl;	
						cout << "repulsive_force of cv_" << i << " = " << rep_force << endl;
						cout << "disp = " << disp_x << ", " << disp_y << endl;
						cout << "node disp = " << node_disp.at(i).first << ", " << node_disp.at(i).second << endl;
					} 
					
				}
			}

		}

		// step 2: attractive force placement
		for (int e=0; e<c_edges.size(); ++e)
		{
			e_u = c_edges.at(e).first;
			e_v = c_edges.at(e).second;

			// \delta <- e.u.pos - e.v.pos
			pos_diff = make_pair(center_coord.at(e_u).at(0)-center_coord.at(e_v).at(0),
								 center_coord.at(e_u).at(1)-center_coord.at(e_v).at(1));

			att_force = attractive_force(pos_diff, radii.at(e_u), radii.at(e_v));

			// e.u.disp <- e.u.disp - unit-direction * att_force(abs(\delta))
			disp_u_x = node_disp.at(e_u).first - pos_diff.first/pos_diff_dist(pos_diff)*att_force;
			disp_u_y = node_disp.at(e_u).second - pos_diff.second/pos_diff_dist(pos_diff)*att_force;
			node_disp.at(e_u) = make_pair(disp_u_x+node_disp.at(e_u).first, disp_u_y+node_disp.at(e_u).second);

			// e.v.disp <- e.v.disp - unit-direction * att_force(abs(\delta))
			disp_v_x = node_disp.at(e_v).first + pos_diff.first/pos_diff_dist(pos_diff)*att_force;
			disp_v_y = node_disp.at(e_v).second + pos_diff.second/pos_diff_dist(pos_diff)*att_force;
			node_disp.at(e_v) = make_pair(disp_v_x+node_disp.at(e_v).first, disp_v_y+node_disp.at(e_v).second);
			if (t==0)
			{
				cout << "pos diff of cv_" << e << " = " << pos_diff.first << ", " << pos_diff.second << endl;	
				cout << "attractive_force of cv_" << e << " = " << att_force << endl;
				cout << "disp = " << disp_v_x << ", " << disp_v_y << endl;
			} 
		}

		// step 3: torque equilibrium
		fill(rotate_angles.begin(), rotate_angles.end(), 0);
		for (int i=0; i<cg.get_num_vtxs(); ++i)
		{
			cluster_nodes.resize(0);
        	for (int c=0; c<clusters.size(); ++c)
            	if (clusters.at(c) == i) cluster_nodes.push_back(c);
			
			calculate_rotate_angle(g, i, cluster_nodes, clusters,
				coord, center_coord, radii, rotate_angles);

		}

		// step 4: shift and rotate the cluster nodes and corresponding nodes
		for (int i=0; i<cg.get_num_vtxs(); ++i)
		{
			cluster_nodes.resize(0);
        	for (int c=0; c<clusters.size(); ++c)
            	if (clusters.at(c) == i) cluster_nodes.push_back(c);

			disp_val = sqrt(pow(node_disp.at(i).first, 2) + pow(node_disp.at(i).second, 2));
			c_x = center_coord.at(i).at(0);
			c_y = center_coord.at(i).at(1);
			c_x += step * node_disp.at(i).first/disp_val;
			c_y += step * node_disp.at(i).second/disp_val;

			for (int j=0; j<cluster_nodes.size(); ++j)
			{
				angle = rotate_angles.at(i);
				// shift nodes
				coord.at(cluster_nodes.at(j)).at(0) -= center_coord.at(i).at(0);
				coord.at(cluster_nodes.at(j)).at(1) -= center_coord.at(i).at(1);

				// rotate nodes
				new_x = coord.at(cluster_nodes.at(j)).at(0);
        		new_y = coord.at(cluster_nodes.at(j)).at(1);
				coord.at(cluster_nodes.at(j)).at(0) = new_x*cos(angle) + new_y*sin(angle);
        		coord.at(cluster_nodes.at(j)).at(1) = -new_x*sin(angle) + new_y*cos(angle);

				// match to new center
				coord.at(cluster_nodes.at(j)).at(0) += c_x;
				coord.at(cluster_nodes.at(j)).at(1) += c_y;

			}

			center_coord.at(i).at(0) = c_x;
			center_coord.at(i).at(1) = c_y;


		}
		

		// update step
		calculate_energy(cluster_dist_mat, center_coord, curr_energy);
		update_step(step, pre_energy, curr_energy);

		cout << "curr_energy=" << curr_energy << endl;
		cout << "step=" << step << endl;

		// check convergence


		// dump information of each iteration
		if (t < 5)
		{
			// coordinates
			fo.open("out/"+string(outfilePath)+"_"+to_string(t)+".coord", fstream::out);
			for (int i=0; i<coord.size(); ++i)
			{
				fo << coord.at(i).at(0) << " " << coord.at(i).at(1) << endl;
			}
			fo.close();

			// center coord
			fo.open("out/"+string(outfilePath)+"_"+to_string(t)+".center", fstream::out);
			for (int i=0; i<center_coord.size(); ++i)
			{
				fo << center_coord.at(i).at(0) << " " << center_coord.at(i).at(0) << endl;
			}
			fo.close();

			// angle
			fo.open("out/"+string(outfilePath)+"_"+to_string(t)+".angle", fstream::out);
			for (int i=0; i<rotate_angles.size(); ++i)
			{
				fo << rotate_angles.at(i) << endl;
			}
			fo.close();
		}
	}

}



/******************************************************************************
 *                            Main Process                                    *
 ******************************************************************************/

int inter_layout(char* outfilePath, Graph::Graph& g,
	DenseMat& distMat,
    std::vector< std::vector<CoordType> >& center_coord,
    std::vector< WgtType >& radii,
    int cls,
    std::vector<int>& clusters,
    std::vector< std::vector<CoordType> >& coord)
{

	// Steps
	// 1. create cluster graph
	// 2. perform force directed method
	//   1) repulsive
	//   2) attractive
	//   3) torque

    using namespace std;


    // Step 1
    int cg_size = center_coord.size();
    Graph cg(cg_size);
    create_cluster_graph(g, clusters, cls, radii, cg);
    cg.print_graph();

    // Step 2
    DenseMat cluster_dist_mat(cg_size, cg_size);
    distance_matrix(cg, cluster_dist_mat);

    // cout << "cluster distance matrix" << endl;
    // cout << cluster_dist_mat << endl;

    
    // create edges vector
    vector< VtxType > nbors;
    vector< pair<int, int> > c_edges;
    for (int i=0; i<cg.get_num_vtxs(); ++i)
    {
    	nbors = cg.adj(i);
    	for (int n=0; n<nbors.size(); ++n)
    	{
    		if (nbors.at(n) > i) c_edges.push_back(make_pair(i, nbors.at(n)));
    	}
    }

    // cout << "cluster edges" << endl;
    // for (int i=0; i<c_edges.size(); ++i)
    // {
    // 	cout << c_edges.at(i).first << " " << c_edges.at(i).second << endl;
    // }

    double init_step = 0.1*cluster_dist_mat.maxCoeff(); // 0.1*2*distMat.maxCoeff()?


    cout << "previous center_coord" << endl;
    for (int i=0; i<center_coord.size(); ++i)
    {
    	cout << center_coord.at(i).at(0) << ", " << center_coord.at(i).at(1) << endl;
    }

    force_direct_with_torque(outfilePath, g, cg, c_edges, cluster_dist_mat, clusters,
    						 center_coord, radii, coord, 100, init_step);

    cout << "after center_coord" << endl;
    for (int i=0; i<center_coord.size(); ++i)
    {
    	cout << center_coord.at(i).at(0) << ", " << center_coord.at(i).at(1) << endl;
    }


    return SUCCESS_INTER_LAYOUT;
}