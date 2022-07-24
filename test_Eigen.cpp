#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <algorithm>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using namespace Eigen;


MatrixXd matr(int a, int b){
    MatrixXd A(a, a);
    A << 1.0, 0.0,
          0.0, 4.0;
    return A; 
}

int main()
{

     /*
     VectorXd p(3);
     p << 1, 2, 3;
     cout << p[0] << p[1] << p[2] << endl;
     p[0] = 100;
     cout << p[0] << p[1] << p[2] << endl;
     */


     /*
     Eigen::MatrixXf A(2, 2);
     A << 1.0, 2.0,
          3.0, 4.0;
     cout << A << endl;    
     */



     /*
     MatrixXd A = matr(2, 2);
     cout << "A:" << A << endl;

     EigenSolver<MatrixXd> es(A);
     MatrixXcd Es = es.eigenvalues();
     MatrixXd Es_real = Es.real();
     //cout << Es_real << endl;

     MatrixXcd vecs = es.eigenvectors();
     MatrixXd vecs_real = vecs.real();
     //cout << "vecs_real" << vecs_real << endl;
     VectorXd gs = vecs_real.row(0);
     cout << "gs" << gs << endl;
     cout << gs.maxCoeff() << endl;
     */    
 
     /*
     Eigen::MatrixXd A(2, 2);
     A << 1.0, 2.0,
          3.0, 4.0;
     cout << A.row(0) << endl;    
     */

     //cout << pow(4., 2.) << endl;


     /*
     Eigen::MatrixXd A(2, 2);
     A.row(0) << 1, 2;
     A.row(1) << 3, 4;
     //A << 1.0, 2.0,
     //     3.0, 4.0;
     //cout << A.row(1) << endl;

     
     Eigen::MatrixXd B(2, 2);
     B.row(0) << 1, 2;
     B.row(1) << 3, 4;

     bool tr = (A == B);
     cout << tr << endl; 
     */
     
     /*
     VectorXd p(3);
     p << 1,2,3;
     //MatrixXd matr;
     cout << Map<MatrixXd>(p.data(), 1, 3) << endl;  
     */


     Eigen::MatrixXd A(2, 2);
     A.row(0) << 1, 2;
     A.row(1) << 3, 4;
  
     Eigen::MatrixXd B(2, 2);
     B.row(0) = A.row(0);
     //B.row(1) << A.row(1);
 
     cout << B << endl;
     //cout << B.cols() << endl;
     //cout << B(0, 1) << endl;
     B = B * 3;
     cout << typeid(B).name() << endl;
     cout << typeid(B.row(1)).name() << endl;
     cout << B.row(1)[1]
 
     //VectorXd p(3);
     //p << 1,2,3;
     //MatrixXd matr;
     //cout << Map<vector<double>>(p.data(), 1, 3) << endl;  
 

}


