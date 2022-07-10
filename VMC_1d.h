# include <iostream>
# include <vector>
# include <eigen3/Eigen/Eigen>

using namespace std;
using namespace Eigen;


class Hubbard_1d
{
private:
    int Lsite;
    int N;
    double t;
    double U;
    MatrixXd Hfree;
    MatrixXd eigvals;
    MatrixXd eigvecs;

public:
    Hubbard_1d(int LL, int NN, double tt, double UU);
      
    MatrixXd get_Hfree(void);
    MatrixXd get_Vstate(void);

    double x_psi0(vector<int>);     
    //MatrixXd all_hop(vector<int>);
    MatrixXd all_hop(vector<int>);
    int count_double_occ(vector<int>);
    double gutzwiller_weight(vector<int>, double);
    double x_psi(vector<int>, double);
    double Eloc(vector<int>, double);

    vector<int> random_init(void);

    vector<int> jump(vector<int>);


    void metropolis(vector<int> init_state, double g, int nthermal, int npoints, int nacc);
};


