# include <iostream>
# include <eigen3/Eigen/Eigen>
# include <string> 
# include <tuple>
# include <vector>
# include <algorithm>
# include "VMC_1d.h"
# include <cmath>
# include <stdlib.h>
#include<numeric>

using namespace std;
using namespace Eigen;

void print1(vector<int> arr)
{
    for (int i = 0; i < arr.size(); i++) 
    { 
        cout << arr[i] << " ";
    }
    cout << endl;
}


void print2(vector<vector<int>> arr)
{
    for (int i = 0; i < arr.size(); i++)
    {
        print1(arr[i]);
    }
    cout << endl;
}

VectorXd vector_to_Vector(vector<int> vec)
{
    VectorXd Vec(vec.size());
    for (int i=0; i < vec.size(); i++)
    {
        Vec[i] = vec[i];
    }
    return Vec;
}

void mod_n(MatrixXd& matr, const int Lsite)
{
    for (int i=0; i<matr.rows(); i++)
    {
        for (int j=0; j<matr.cols(); j++)
        {
            matr(i, j) = (int) matr(i, j) % Lsite;
            if (matr(i, j) < 0)
            {
                matr(i, j) += Lsite;
            }
        }
    }    
}


vector<int> argsort(vector<double>& arr){
    const int N = arr.size(); 
    vector<int> idx;
    for (int k = 0; k < N; k++)
    {
        idx.push_back(k);
    }

    for (int i=N-1; i>=0; i--)
    {
        for (int j=0; j<=i-1; j++)
        {
            if (arr[j] > arr[j+1])
            {
                double tmp2 = idx[j];
                idx[j] = idx[j+1];
                idx[j+1] = tmp2;
            }
        }
    }   
    return idx;
}



MatrixXd vstack(MatrixXd ma, MatrixXd mb)
{
    int out_rows = ma.rows() + mb.rows();
    int out_cols = ma.cols();
    MatrixXd out(out_rows, out_cols);
    for (int i=0; i < ma.rows(); i++)
    {
        out.row(i) = ma.row(i);
    }
    for (int j=0; j < mb.rows(); j++)
    {
        out.row(ma.rows() +j) = mb.row(j);
    }
    return out;
}


MatrixXd hstack(MatrixXd ma, MatrixXd mb)
{
    int out_rows = ma.rows();
    int out_cols = ma.cols() + mb.cols();
    MatrixXd out(out_rows, out_cols);
    for (int i=0; i < ma.cols(); i++)
    {
        out.col(i) = ma.col(i);
    }
    for (int j=0; j < mb.cols(); j++)
    {
        out.col(ma.cols() +j) = mb.col(j);
    }
    return out;
}


/*
MatrixXd hstack(MatrixXd ma, MatrixXd mb)
{
    int out_rows = ma.rows();
    int out_cols = ma.cols() + mb.cols();
    MatrixXd out(out_rows, out_cols);
    for (int i=0; i<out_rows; i++)
    {
        for (int j=0; j < ma.cols(); j++)
        {
            out(i, j) = ma(i, j);
        }
        for (int k=0; k < mb.cols(); k++)
        {
            out(i, k+ma.cols()) = mb(i, k);
        }
    }
    return out;
}
*/

// class memember functions ========================================================================


Hubbard_1d::Hubbard_1d(int LL, int NN, double tt, double UU){
    Lsite = LL;
    N = NN;
    t = tt;
    U = UU;

        MatrixXd H_free(Lsite*2, Lsite*2);
        for (int i=1; i < Lsite-1; i++)
        {
            H_free(i, i+1) = -t;
            H_free(i, i-1) = -t;
        }
        H_free(0, 1) = -t;
        H_free(0, Lsite-1) = -t;
        H_free(Lsite-1, Lsite-2) = -t;
        H_free(Lsite-1, 0) = -t;

        for (int j=Lsite+1; j < Lsite*2-1; j++)
        {
            H_free(j, j+1) = -t;
            H_free(j, j-1) = -t;
        }
        H_free(Lsite, Lsite+1) = -t;
        H_free(Lsite, Lsite*2-1) = -t;
        H_free(Lsite*2-1, Lsite*2-2) = -t;
        H_free(Lsite*2-1, Lsite) = -t;
  
        Hfree = H_free;

        SelfAdjointEigenSolver<MatrixXd> es(H_free);
        MatrixXd Es = es.eigenvalues();
        MatrixXd vecs = es.eigenvectors();
         
        eigvals = Es;
        eigvecs = vecs;
}             
 

MatrixXd Hubbard_1d::get_Vstate(void)
{
    MatrixXd Vstate(Lsite*2, N*2); 
    for (int i=0; i < N*2; i++)
    {
        Vstate.col(i) = eigvecs.col(i);
    }
    return Vstate;
}

double Hubbard_1d::x_psi0(vector<int> x)
{
    MatrixXd VV = get_Vstate();
    MatrixXd VVx(x.size(), N*2);
    for (int i=0; i < N*2; i++)
    {
        VVx.row(i) = VV.row(x[i]);
    }
    return VVx.determinant();
}

MatrixXd Hubbard_1d::all_hop(vector<int> x)
{
    vector<int> up, down;
    for (int i=0; i < N; i++)
    {
        up.push_back(x[i]);
        down.push_back(x[i+N] - Lsite);
    }
    MatrixXd up_copy(N*2, N), down_copy(N*2, N);
  
    for (int i=0; i < N*2; i++)
    {
        up_copy.row(i) = vector_to_Vector(up);
        down_copy.row(i) = vector_to_Vector(down);
    } 

    MatrixXd iden = MatrixXd::Identity(N, N);
    MatrixXd matr_hop = vstack(iden, -iden);

    MatrixXd up_hopped = up_copy + matr_hop; 
    MatrixXd down_hopped = down_copy + matr_hop; 

    mod_n(up_hopped, Lsite);
    mod_n(down_hopped, Lsite);

    MatrixXd ones = MatrixXd::Ones(N*2, N);
    down_copy += Lsite * ones;
    down_hopped += Lsite * ones;

    MatrixXd hopped_up = hstack(up_hopped, down_copy);
    MatrixXd hopped_down = hstack(up_copy, down_hopped);

    return vstack(hopped_up, hopped_down);
}

int Hubbard_1d::count_double_occ(vector<int> state)
{
    int double_occ = 0;
    for (int i=0; i<N; i++)
    {
        int key = (state[i] + Lsite);
        if (std::count(state.begin()+N, state.end(), key))
        {
            double_occ += 1;
        }
    }     
    return double_occ;
}

double Hubbard_1d::gutzwiller_weight(vector<int> state, double g)
{
    int double_occ = count_double_occ(state);    
    return pow(g, double_occ);
}

double Hubbard_1d::x_psi(vector<int> x, double g)
{
    return x_psi0(x) * gutzwiller_weight(x, g);
}

double Hubbard_1d::Eloc(vector<int> x, double g)
{
    MatrixXd all_hopped = all_hop(x);
    //cout << all_hopped << endl;
    
    double kinetic = 0.;
     
    for (int i=0; i < all_hopped.rows(); i++)
    {

        //if (i != 3)
        //{
        //    continue;
        //}

        vector<int> x_hopped;
        for (int j=0; j < all_hopped.cols(); j++)
        {
            x_hopped.push_back(all_hopped(i, j));
        }
        //print1(x_hopped);

        double gutzwiller_ratio = gutzwiller_weight(x_hopped, g) / gutzwiller_weight(x, g);
        //cout << "gutz_ratio:" << " " << gutzwiller_ratio << endl;

        double psi0_ratio = x_psi0(x_hopped) / x_psi0(x);
        //cout << x_psi0(x_hopped) << endl;

        double psi_ratio = psi0_ratio * gutzwiller_ratio;
        //cout << "psi_ratio:" << " " << psi_ratio << endl;
        kinetic += (-t) * psi_ratio;
    }
    double potential = U * count_double_occ(x);
    return kinetic + potential;
}


vector<int> Hubbard_1d::random_init(void)
{ 
    vector<int> site_all;
    for (int i=0; i < Lsite; i++)
    {
        site_all.push_back(i);
    } 
    std::random_shuffle(site_all.begin(), site_all.end());
    //print1(site_all);
 
    vector<int> site_up;
    for (int i=0; i < N; i++)
    {
        site_up.push_back(site_all[i]);
    } 

    std::random_shuffle(site_all.begin(), site_all.end());
    print1(site_all);
 
    vector<int> site_down;
    for (int i=0; i < N; i++)
    {
        site_down.push_back(site_all[i]+Lsite);
    }

    vector<int> state_init;
    for (int i=0; i < N; i++)
    {
        state_init.push_back(site_up[i]);
    }
    for (int j=0; j < N; j++)
    {
        state_init.push_back(site_down[j]);
    }
    return state_init;
}


vector<int> Hubbard_1d::jump(vector<int> state)
{
    int key_spin = rand() % 2;
    int key_x = rand() % N;
    int key_proposal = rand() % N;

    //cout << "key_spin:" << " " << key_spin << endl;
    //cout << "key_x:" << " " << key_x << endl;
    //cout << "key_proposal:" << " " << key_proposal << endl;

    vector<int> x_proposal;

    if (key_spin == 0) // spin up
    {
        for (int i=0; i < Lsite; i++)
        {
            bool tr = std::count(state.begin(), state.end(), i);
            if (tr==0)
            {
                x_proposal.push_back(i);
            }
        }
        state[key_x] = x_proposal[key_proposal];
    }
    else  // spin down
    {
        for (int j=Lsite; j < Lsite*2; j++)
        {
            bool tr = std::count(state.begin(), state.end(), j);
            if (tr==0)
            {
                x_proposal.push_back(j);
            }
        }
        state[N+key_x] = x_proposal[key_proposal];
    }
    //cout << "x_proposal:" << endl;
    //print1(x_proposal);
    return state; 
}


void Hubbard_1d::metropolis(vector<int> init_state, double g, int nthermal, int npoints, int nacc)
{
    int nstep = nthermal + npoints * nacc;
    vector<double> eloc_all;

    vector<int> state = init_state;

    int ipoint = 0;
    for (int istep = 0; istep < nstep; istep++)
    {
        vector<int> state_proposal = jump(state);
        double ratio = x_psi(state_proposal, g) / x_psi(state, g);
        double accept_rate = pow(ratio, 2);
        double eta = rand() / double(RAND_MAX);
        cout << "eta:" << " " << eta << endl;
        
        if (eta <= accept_rate)
        {
            state = state_proposal;
        }        
        double eloc = Eloc(state, g);
        cout << "eloc:" << " " << eloc << endl;
        eloc_all.push_back(eloc);
    }    

    double E_sum = std::accumulate( eloc_all.begin(), eloc_all.end(), 0.0 );
    double E_mean = E_sum / eloc_all.size();
    cout << "E_mean:" << " " << E_mean << endl;
    cout << "E_mean_per_site:" << " " << E_mean/Lsite << endl;
}


// test ========================================================================

void test_hstack()
{
    MatrixXd A(3, 5), B(3, 2);
    cout << hstack(A, B) << endl;
}


/*
void test_Hfree()
{
    Hubbard_1d model = Hubbard_1d(6, 3, 1., 6.);
    cout << model.get_Hfree() << endl;
    cout << endl;
}
*/

void test_Vstate()
{
    //Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    Hubbard_1d model = Hubbard_1d(6, 3, 1., 6.);
    for (int i=0; i < 10; i++)
    {
        cout << model.get_Vstate() << endl;
        cout << endl;
    }
}

void test_x_psi0()
{
    //Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    //vector<int> x = {0, 1, 4, 5};
    Hubbard_1d model = Hubbard_1d(6, 3, 1., 6.);
    vector<int> x = {5, 1, 2, 6, 7, 8};
    for (int i=0; i < 10; i++)
    {
        cout << model.x_psi0(x) << endl;
        cout << endl;
    }
}

void test_all_hop()
{
    Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    vector<int> x = {0, 1, 4, 5};
    model.all_hop(x);
}

void test_gutzwiller_weight()
{
    Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    vector<int> x = {0, 1, 4, 5};
    double g = 0.5;
    cout << model.gutzwiller_weight(x, g) << endl;
}

void test_Eloc()
{
    //Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    //vector<int> x = {1, 2, 4, 5};
    Hubbard_1d model = Hubbard_1d(6, 3, 1., 6.);
    vector<int> x = {0, 1, 2, 6, 7, 8};
    double g = 0.5;
    for (int i=0; i < 50; i++)
    {
        cout << model.Eloc(x, g) << endl;
        cout << endl;
    }
}


void test_random_init()
{
    //Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    //vector<int> x = {1, 2, 4, 5};
    Hubbard_1d model = Hubbard_1d(6, 3, 1., 6.);
    for (int i=0; i < 50; i++)
    {
        print1( model.random_init() );
        cout << endl;
    }
}


void test_jump()
{
    Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    vector<int> x = {0, 1, 4, 5};
    //Hubbard_1d model = Hubbard_1d(6, 3, 1., 6.);
    for (int i=0; i < 10; i++)
    {
        //print1( model.jump(x) );
        vector<int> x_new = model.jump(x);
        print1(x_new);
        print1(x);
        cout << endl;
    }
}


void test_metropolis()
{
    //Hubbard_1d model = Hubbard_1d(4, 2, 1., 6.);
    //vector<int> x = {0, 1, 4, 5};
    Hubbard_1d model = Hubbard_1d(30, 15, 1., 4.);
    //vector<int> x = {0, 1, 4, 5};
    vector<int> state = model.random_init();
    double g = 0.47;
    int nthermal = 200;
    int npoints = 1000;
    int nacc = 1;
    model.metropolis(state, g, nthermal, npoints, nacc);
}

// ========================================================================


int main()
{
    //test_hstack();
    //vector<int> v = {1, 2, 3};
    //print1(v);
    //cout << vector_to_Vector(v) << endl;
    //vector<vector<int>> vv = {{1, 2, 3}, {4, 5, 6}};
    //print2(vv);

    //test_Hfree();
    //test_Vstate();
    //test_x_psi0();

    /*
    for (int i = 0; i < 10; i ++)
    {
        cout << i << endl;
        //test_Hfree();
        test_Vstate();
        //test_x_psi0();
        //test_Eloc();
    }
    */


    //test_all_hop();
    //test_gutzwiller_weight();
    //test_Eloc();
    //test_random_init();
    //test_jump();
    test_metropolis();
}

