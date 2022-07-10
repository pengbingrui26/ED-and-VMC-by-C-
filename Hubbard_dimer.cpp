#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <cmath>
using namespace std;
using Eigen::MatrixXd;
using namespace Eigen;


class Dimer {
public:
    double t;
    double U; 
    MatrixXd H;

    Dimer(double tt, double UU){
    t = tt;
    U = UU;
    MatrixXd H_tmp(4, 4);
    H_tmp << 0., 0., -t, -t,
         0., 0., -t, -t,
         -t, -t, U, 0.,
         -t, -t, 0., U;
    H = H_tmp;
    }

    MatrixXd get_eigvals(void);
    MatrixXd get_eigvecs(void);
    void show_eigs(void);
    VectorXd get_GS(void);
    MatrixXd get_Hfree(void);
    VectorXd get_gwf(double g);
    double partition_fn(double beta);    
    double free_E(double beta);    
};


MatrixXd Dimer::get_Hfree(void){
    Dimer model_tmp = Dimer(t, 0.);
    return model_tmp.H;
}

void Dimer::show_eigs(void){
    SelfAdjointEigenSolver<MatrixXd> es(H);
    MatrixXd Es = es.eigenvalues();
    MatrixXd vecs = es.eigenvectors();
    cout << Es << endl; 
    cout << vecs << endl; 
}


//Column k of the returned matrix is an eigenvector corresponding to eigenvalue number k as returned by eigenvalues()
MatrixXd Dimer::get_eigvals(void){
    VectorXd E(4);
    SelfAdjointEigenSolver<MatrixXd> es(H);
    MatrixXd Es = es.eigenvalues();
    return Es;
}

MatrixXd Dimer::get_eigvecs(void){
    SelfAdjointEigenSolver<MatrixXd> es(H);
    MatrixXd vecs = es.eigenvectors();
    return vecs.real();
}

//get ground state
VectorXd Dimer::get_GS(void){
    MatrixXd vecs = get_eigvecs();
    return vecs.col(0);
}


VectorXd Dimer::get_gwf(double g){
    MatrixXd H_free = get_Hfree();
    SelfAdjointEigenSolver<MatrixXd> es(H_free);
    VectorXd gs_free = es.eigenvectors().col(0);
    //cout << "gs_free" << gs_free << endl;
    MatrixXd double_occ(4,4);    
    double_occ << pow(g, 0), 0, 0, 0,
         0, pow(g, 0), 0, 0,
         0, 0, pow(g, 1), 0,
         0, 0, 0, pow(g, 1);
    //cout << double_occ << endl;
    VectorXd gwf = double_occ * gs_free;
    //cout << gwf << endl;
    return gwf;
}

double Dimer::partition_fn(double beta){
    MatrixXd Es = get_eigvals();
    cout << Es << endl;
    const int N = Es.size();
    double Z = 0.;
    for (int i=0; i<N; i++)
    {
        Z += exp(-beta * Es(i,0));
    } 
    return Z;
}


double Dimer::free_E(double beta){
    MatrixXd Es = get_eigvals();
    cout << Es << endl;
    const int N = Es.size();
    double F = 0.;
    double Z = partition_fn(beta);
    for (int i=0; i<N; i++)
    {
        double p = exp(-beta * Es(i,0))/Z;
        double f = 1/beta * p * log(p) + p * Es(i,0);
        F += f;
    } 
    return F;
}


int main(){
    Dimer model = Dimer(1, 6.); 
    //cout << model.t << endl;
    //cout << "model.H" << model.H << endl;
    //cout << "model.H_free" << model.get_Hfree() << endl;

    //model.show_eigs();

    MatrixXd Es = model.get_eigvals();
    //cout << Es << endl;
    MatrixXd vecs = model.get_eigvecs();
    //cout << vecs << endl;
    VectorXd gs = model.get_GS();
    //cout << "gs:" << gs << endl;

    //cout << model.get_gwf(0.5) << endl;
    cout << model.free_E(1.) << endl;
}

