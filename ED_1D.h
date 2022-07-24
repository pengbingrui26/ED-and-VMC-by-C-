// ====================================================================
// Exact diagonalization of 1D Hubabbard model
// ====================================================================

# include<iostream>
# include <eigen3/Eigen/Eigen>

using namespace std;
using namespace Eigen;

class Hubbard_1d{
public:
    int Lsite;
    int N;
    double t;
    double U;
    Hubbard_1d(int LL, int NN, double tt, double UU){
        Lsite = LL;
        N = NN;
        t = tt;
        U = UU;
    }    
    vector<vector<int>> get_basis(void);
    int count_double_occ(vector<int>);
    tuple<vector<int>, int> hop(vector<int>, int x, int dire);

    tuple< vector<vector<int>>, vector<int> > all_hop(vector<int>);
    //void all_hop(vector<int>);

    MatrixXd get_T(void);
    MatrixXd get_U(void);
    MatrixXd get_H(void);
    MatrixXd get_eigvals(void);
    MatrixXd get_eigvecs(void);

};




