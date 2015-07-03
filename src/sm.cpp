#include "sm.h"
#include <ctime>


static
double compute_norm(DenseMat& coord, int u, int v)
{
    using namespace std;
    double norm = 0;
    norm += (coord(u, 0) - coord(v, 0)) * (coord(u, 0) - coord(v, 0));
    norm += (coord(u, 1) - coord(v, 1)) * (coord(u, 1) - coord(v, 1));
    return sqrt(norm);
}


static
double stress(DenseMat & distMat, DenseMat& coord)
{
    using namespace std;
    double stress = 0.0;
    int nNodes = coord.rows();
    double dist;
    double diff;

    // other node not central node
    for (int i=0; i<nNodes; ++i)
    {
        for (int j=i+1; j<nNodes; ++j)
        {
            if (i != j)
            {
                dist = distMat(i, j);
                diff = compute_norm(coord, i, j) - dist;
                stress += pow(dist, -2) * pow( diff, 2 );
            }
        }
    }

    return stress;
}


static
double stress_with_radial(DenseMat & distMat, DenseMat& coord, double coeff)
{
    using namespace std;
    double stress = 0.0;
    int nNodes = coord.rows();
    // radial part
    for (int i=1; i<nNodes; ++i)
    {
        stress += coeff * pow(distMat(0, i), -2) * \
            pow( ((coord.row(0) - coord.row(i)).norm() - distMat(0, i)), 2 );
    }
    // other node not central node
    for (int i=1; i<nNodes; ++i)
    {
        for (int j=i+1; j<nNodes; ++j)
        {
            if (i != j)
            {
                stress += (1-coeff) * pow(distMat(i, j), -2) * \
                pow( ((coord.row(i) - coord.row(j)).norm() - distMat(i, j)), 2 );
            }
        }
    }

    return stress;
}


int stress_majorization(int graphSize,
                        DenseMat& distMat,
                        std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Create Corresponding Coordinates Matrix
    // 2. Create weight laplacian matrix
    // 3. Majorization
    //     * Pick the first node as reference point
    // 4. map Eigen matrix to coord
    using namespace std;
    clock_t t_sm_begin = clock();


    int g_size = graphSize;

    // Step 1
    DenseMat coord_mat(g_size, 2);
    for (int c=0; c<coord_mat.cols(); c++)
    {
        for (int r=0; r<coord_mat.rows(); ++r)
        {
            coord_mat(r, c) = coord.at(r).at(c);
        }
    }
    

    // Step 2
    DenseMat w_lap(g_size, g_size);
    w_lap_normal(g_size, distMat, w_lap);


    // Step 3
    double stress_ratio = 0.0;
    double epsl = SM_THRESHOLD;
    double pre_stress;
    double aft_stress;
    DenseMat i_lap(w_lap.rows(), w_lap.cols());
    DenseMat p_sol;
    DenseMat iter_coord;
    double dist_norm;

    Eigen::ConjugateGradient< DenseMat > conj_g; // for conjugate gradient

    int r_cnt = 0;
    while (true)
    {
        clock_t t_stress_cal_begin = clock();

        pre_stress = stress(distMat, coord_mat);

        clock_t t_stress_cal_end = clock();
        double t_stress_cal = double(t_stress_cal_end - t_stress_cal_begin) / CLOCKS_PER_SEC;
        if (r_cnt == 0) cout << "stress calculation elapsed time = " << t_stress_cal << " secs" << endl;
        cout << "previous stress=" << pre_stress << endl;

        // create i_lap
        clock_t t_i_lap_fill_begin = clock();

        i_lap.fill(0);


        clock_t t_i_lap_fill_end = clock();
        double t_i_lap_fill = double(t_i_lap_fill_end - t_i_lap_fill_begin) / CLOCKS_PER_SEC;
        if (r_cnt == 0) cout << "i lap fill elapsed time = " << t_i_lap_fill << " secs" << endl;


        clock_t t_i_lap_create_begin = clock();

        for (int r=0; r<i_lap.rows(); ++r)
        {
            
            for (int c=0; c<i_lap.cols(); ++c)
            {

                if (r != c)
                {
                    // dist_norm = (coord_mat.row(r) - coord_mat.row(c)).norm();
                    dist_norm = compute_norm(coord_mat, r, c);

                    if (dist_norm==0) dist_norm = 1; // [Warning] incase norm = 0;
                    i_lap(r, c) = w_lap(r, c) * distMat(r, c) / dist_norm;
                    i_lap(r, r) -= i_lap(r, c);
                }
            }
        }

        clock_t t_i_lap_create_end = clock();
        double t_i_lap_create = double(t_i_lap_create_end - t_i_lap_create_begin) / CLOCKS_PER_SEC;
        if (r_cnt == 0) cout << "i lap create elapsed time = " << t_i_lap_create << " secs" << endl;


        // solving part

        clock_t t_p_sol_multi_begin = clock();

        p_sol = i_lap * coord_mat;

        clock_t t_p_sol_multi_end = clock();
        double t_p_sol_multi = double(t_p_sol_multi_end - t_p_sol_multi_begin) / CLOCKS_PER_SEC;
        if (r_cnt == 0) cout << "p sol multiplication elapsed time = " << t_p_sol_multi << " secs" << endl;


        clock_t t_iter_solve_begin = clock();

        // // Cholesky Decomposition
        // iter_coord = w_lap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1)\
        //     .ldlt().solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));

        // Conjugate Gradient
        conj_g.compute(w_lap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1));
        iter_coord = conj_g.solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));


        clock_t t_iter_solve_end = clock();
        double t_iter_solve = double(t_iter_solve_end - t_iter_solve_begin) / CLOCKS_PER_SEC;
        if (r_cnt == 0)  cout << "iter solve elapsed time = " << t_iter_solve << " secs" << endl;


        coord_mat.fill(0);  // refill coord_mat for next iteration


        clock_t t_coord_block_begin = clock();

        coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;

        clock_t t_coord_block_end = clock();
        double t_coord_block = double(t_coord_block_end - t_coord_block_begin) / CLOCKS_PER_SEC;
        if (r_cnt == 0) cout << "coord block elapsed time = " << t_coord_block << " secs" << endl;


        aft_stress = stress(distMat, coord_mat);
        // cout << "after stress=" << aft_stress << endl;
        // decide whether to stop stress majorization
        stress_ratio = (pre_stress-aft_stress)/pre_stress;

        r_cnt++;
        if ( (stress_ratio < epsl) || isnan(stress_ratio) ) break;

    }
    cout << "round count = " << r_cnt << endl;


    // transform coord_mat back to pg_coord
    for (int c=0; c<coord_mat.cols(); c++)
    {
        for (int r=0; r<coord_mat.rows(); ++r)
        {
        
            coord.at(r).at(c) = coord_mat(r, c);
        }
    }
    clock_t t_sm_end = clock();
    double t_sm = double(t_sm_end - t_sm_begin) / CLOCKS_PER_SEC;
    cout << "stress_majorization elapsed time = " << t_sm << " secs" << endl;

    return SUCCESS_SM;
}


int stress_majorization_radial_refinement(
    DenseMat& distMat,
    double coeff,   // the coeff of linear combs of orig stress and constriants
    std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Create weight laplacian matrix
    // 2. Majorization
    //     * Pick the first node (center) as reference point
    // 3. map Eigen matrix to coord
    using namespace std;
    int coord_size = coord.size();

    
    DenseMat coord_mat(coord_size, 2);
    int coord_mat_cols = coord_mat.cols();
    int coord_mat_rows = coord_mat.rows();
    for (int c=0; c<coord_mat_cols; c++)
    {
        for (int r=0; r<coord_mat_rows; ++r)
        {
            coord_mat(r, c) = coord.at(r).at(c);
        }
    }
    // cout << "init coord" << endl;
    // cout << coord_mat << endl;



    // Step 1
    DenseMat w_lap(coord_size, coord_size);
    w_lap.fill(0);
    int w_lap_size = w_lap.rows();

    // deal the radius part
    for (int i=1; i<w_lap_size; ++i) {
        w_lap(0, i) = -1 * coeff * pow(distMat(0, i), -2);
        w_lap(i, 0) = -1 * coeff * pow(distMat(i, 0), -2);

        // i == j part
        w_lap(i, i) -= w_lap(0, i);
    }

    // deal the i != j other part
    for (int i=1; i<w_lap_size; ++i) {

        // i != j
        for (int j=1; j<w_lap_size; ++j)
        {
            if (i != j)
            {
                w_lap(i, j) = -1 * (1-coeff) * pow(distMat(i, j), -2);
                w_lap(i, i) -= w_lap(i, j); // i == j part
            } 

        }

    }


    // cout << "distance matrix" << endl;
    // cout << distMat << endl;
    // cout << "w lap" << endl;
    // cout << w_lap << endl;


    // Step 3
    double stress_ratio = 0.0;
    double dist_norm;
    double epsl = SM_THRESHOLD;
    double pre_stress;
    double aft_stress;
    DenseMat i_lap(w_lap.rows(), w_lap.cols());
    DenseMat p_sol;
    DenseMat iter_coord;
    Eigen::ConjugateGradient< DenseMat > conj_g; // for conjugate gradient

    cout << "coeff=" << coeff << endl;

    int i = 0;
    while (true)
    {
        pre_stress = stress_with_radial(distMat, coord_mat, coeff);
        cout << "previous stress=" << pre_stress << endl;
        // create i_lap
        i_lap.fill(0);
        for (int r=0; r<i_lap.rows(); ++r)
        {
            for (int c=0; c<i_lap.cols(); ++c)
            {
                if (r != c)
                {
                    dist_norm = compute_norm(coord_mat, r, c);
                    if (dist_norm==0) dist_norm = 1; // [Warning] incase norm = 0;
                    i_lap(r, c) = w_lap(r, c) * distMat(r, c) / dist_norm;
                    i_lap(r, r) -= i_lap(r, c);
                }
            }
        }

        // solving part
        // cout << "i_lap=" << i_lap << endl;
        p_sol = i_lap * coord_mat;
        // cout << "psol=" << p_sol << endl;

        iter_coord = w_lap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1)\
            .ldlt().solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));
        // cout << "iter_coord=" << iter_coord << endl;

        // conj_g.compute(w_lap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1));
        // iter_coord = conj_g.solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));

        coord_mat.fill(0);  // refill coord_mat for next iteration
        coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
        aft_stress = stress_with_radial(distMat, coord_mat, coeff);
        stress_ratio = (pre_stress-aft_stress)/pre_stress;
        if ( (stress_ratio < epsl) || isnan(stress_ratio) ) break;
        ++i;
    }

    // Step 3
    for (int c=0; c<coord_mat.cols(); c++)
    {
        for (int r=0; r<coord_mat.rows(); ++r)
        {
        
            coord.at(r).at(c) = coord_mat(r, c);
        }
    }


    return SUCCESS_SM;
}