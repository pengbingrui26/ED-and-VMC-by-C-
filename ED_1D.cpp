# include <iostream>
# include <eigen3/Eigen/Eigen>
# include <string> 
# include <tuple>
# include <algorithm>
# include "ED_1D.h"

using namespace std;
using namespace Eigen;


VectorXd vector_to_Vector(vector<int> vec){
    int N = vec.size();
    VectorXd Vec(N);
    for (int i=0; i < vec.size(); i++)
    {
        Vec[i] = vec[i];  
    }
    return Vec;
}


MatrixXd Vector_to_Matrix(VectorXd Vec){
    int N = Vec.size();
    MatrixXd matr = Map<MatrixXd>(Vec.data(), 1, N);
    return matr;
}


void print1(vector<int> arr)
{
    for (int i = 0; i < arr.size(); i++) 
    { 
        cout << arr[i];
    }
    cout << endl;
    //cout << Vector_to_Matrix(vector_to_Vector(arr)) << endl;
}


void comb(vector<vector<int>>& empty, vector<int>& arr, const int n, int m, const int M)
{
    if (m<=1) {
        for (int i=0; i <= n-M; i++){
            vector<int> tmp = {arr[i]};
            empty.push_back(tmp);
        }
    }
    else
    {
        comb(empty, arr, n, m-1, M);
        vector<vector<int>> new_arr;
        for (int j=0; j < empty.size(); j++)
        {
            int lastnum = empty[j].back();
            for (int k=(lastnum+1); k <= arr.back(); k++)
            {
                vector<int> tmp_arr = empty[j];
                tmp_arr.push_back(k);
                new_arr.push_back(tmp_arr);
            }
        }
        empty = new_arr;
    }
}


vector<vector<int>> get_comb(vector<int>& arr, const int M){
    const int n = arr.size();
    int m = M;
    vector<vector<int>> empty;
    comb(empty, arr, n, m, M);
    return empty;
}


void print3(vector<vector<int>> arr)
{
    for (int i = 0; i < arr.size(); i++)
    {
        //cout << vector_to_Vector(arr[i]) << endl; 
        //cout << "\n" << endl;
        print1(arr[i]);
    }
}

void print4(vector<vector<int>> arr)
{
    for (int i = 0; i < arr.size(); i++)
    {
        //cout << Vector_to_Matrix(vector_to_Vector(arr[i])) << endl; 
        print1(arr[i]);
    }
}


vector<int> combine_two(vector<int> arr1, vector<int> arr2){
    vector<int> vec;
    for (int i=0; i<arr1.size(); i++){
        vec.push_back(arr1[i]);
    }
    for (int j=0; j<arr2.size(); j++){
        vec.push_back(arr2[j]);
    }
    return vec;
}


int bubble_sort_parity(vector<int>& arr){
    const int N = arr.size(); 
    int parity = 0;
    for (int i=N-1; i>=0; i--)
    {
        for (int j=0; j<=i-1; j++)
        {
            if (arr[j] > arr[j+1])
            {
                double tmp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp;  
                parity += 1;
            }
        }
    }   
    return parity;
} 

bool vector_in_vvector(vector<int> arr, vector<vector<int>> barr){
    bool tr = 0;
    for (int i = 0; i < barr.size(); i++)
    {
        if (arr == barr[i])
        {
            tr = 1;
            break;
        }
    }
    return tr;
}


// member functions =============================================================

vector<vector<int>> Hubbard_1d::get_basis(void){
    vector<int> sites;
    for (int l=0; l < Lsite; l++)
    {
        sites.push_back(l);
    }
    vector<vector<int>> basis_up = get_comb(sites, N) ;
    vector<vector<int>> basis_down = basis_up;

    vector<vector<int>> basis;
    for (int i=0; i < basis_up.size(); i++)
    {
        for (int j=0; j < basis_down.size(); j++)
        {
             vector<int> up = basis_up[i];
             vector<int> down = basis_down[j];
             for (int s = 0; s < down.size(); s++)
             {
                 down[s] += Lsite;
             }
             vector<int> state = combine_two(up, down);
             //print1(state);
             basis.push_back(state);
        }
    }   
    return basis;
}

int Hubbard_1d::count_double_occ(vector<int> state){
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


tuple<vector<int>, int> Hubbard_1d::hop(vector<int> state, int x, int dire){
    vector<int> state_new;
    int parity = -1;
    string spin;
    if (state[0] < Lsite ){ 
        spin = "up"; 
    } 
    else{ spin = "down";
    }
    if (spin == "down")
    {
        for (int i = 0; i < state.size(); i++)
        {
            state[i] -= Lsite;
        }
    }
    int x11 = (x+dire) % Lsite;
    if (x11 < 0)
    {
        x11 += Lsite;
    }
    bool tr1 = std::count(state.begin(), state.end(), x);
    bool tr2 = std::count(state.begin(), state.end(), x11);
    if (tr1 && ! tr2) 
    {
        //cout << "go on" << endl;
        //print1(state);
        state_new = state;
        for (int k = 0; k < state_new.size(); k++)
        {
            if (state_new[k] == x)
            { 
                state_new[k] = x11;
            }
        }
        parity = bubble_sort_parity(state_new) % 2; 
        if (spin == "down")
        {
            for (int j = 0; j < state_new.size(); j++)
            {
                state_new[j] += Lsite;
            }
        }

    }
    tuple<vector<int>, int> tu(state_new, parity);
    return tu;
}


tuple< vector<vector<int>>, vector<int> > Hubbard_1d::all_hop(vector<int> state){
    vector<int> up;
    vector<int> down;
    for (vector<int>::iterator it = state.begin(); it != state.end()-N; it++)
    {
        //cout << *it << endl;        
        up.push_back(*it);
        down.push_back(*(it + N));        
    }
 
    vector<vector<int>> state_hopped;
    vector<int> parities;
    
    for (int i=0; i < N; i++)
    {
        int x_up = up[i];
        int x_down = down[i] - Lsite;
        for (int dir=-1; dir <=1; dir=dir+2)
        {
            //cout << dir << endl;
            tuple<vector<int>, int> tu_up = hop(up, x_up, dir);
            tuple<vector<int>, int> tu_down = hop(down, x_down, dir);
            vector<int> up_hopped = combine_two( std::get<0>(tu_up), down );
            vector<int> down_hopped = combine_two( up, std::get<0>(tu_down) );

            int parity_up = std::get<1>(tu_up);
            int parity_down = std::get<1>(tu_down);

            //print4(state_hopped);
            //print1(up_hopped);
            //print1(down_hopped);
            //cout << "\n" << endl;

            bool tr1 = vector_in_vvector(up_hopped, state_hopped);
            bool tr2 = vector_in_vvector(down_hopped, state_hopped);

            if (tr1==0 && up_hopped.size() == N*2 )
            {
                state_hopped.push_back(up_hopped);
                parities.push_back(parity_up); 
            }
            if (tr2==0 && down_hopped.size() == N*2 )
            {
                state_hopped.push_back(down_hopped);
                parities.push_back(parity_down); 
            }
 
        }
    //print4(state_hopped);
    //print1(parities);
    } 
    tuple< vector<vector<int>>, vector<int> > tu(state_hopped, parities);
    return tu;
}


MatrixXd Hubbard_1d::get_T(void){
    vector<vector<int>> basis = get_basis();
    //print4(basis);
    //cout << "end of basis" << endl;

    MatrixXd T_matr(basis.size(), basis.size());

    for (int i = 0; i < basis.size(); i++)
    {
        vector<int> state = basis[i];
        tuple< vector<vector<int>>, vector<int> > tu;
        tu = all_hop(state);

        vector<vector<int>> state_hopped = std::get<0>(tu);
        vector<int> parities = std::get<1>(tu);
        //print4(state_hopped);        
 
        for (int j=0; j < state_hopped.size(); j++)
        {
             vector<int> state_new = state_hopped[j];
             int parity = parities[j];
             
             vector<vector<int>>::iterator ip = find(basis.begin(), basis.end(), state_new);
             int idx = std::distance(basis.begin(), ip);
             //print1(state_new);
             //cout << idx << endl;
             //cout << parity << endl;
             T_matr(i, idx) = -t * pow(-1, parity);
 
        }
    }
    //cout << T_matr << endl;
    return T_matr;
}   


MatrixXd Hubbard_1d::get_U(void){
    vector<vector<int>> basis = get_basis();
    //print4(basis);
    //cout << "end of basis" << endl;

    MatrixXd U_matr(basis.size(), basis.size());

    for (int i = 0; i < basis.size(); i++)
    {
        vector<int> state = basis[i];
        int double_occ = count_double_occ(state);
        U_matr(i, i) = U * double_occ;
    }
    return U_matr;
}

MatrixXd Hubbard_1d::get_H(void){
    return get_T() + get_U();
}


MatrixXd Hubbard_1d::get_eigvals(void){
    MatrixXd H = get_H();
    SelfAdjointEigenSolver<MatrixXd> es(H);
    MatrixXd Es = es.eigenvalues();
    return Es;
}

MatrixXd Hubbard_1d::get_eigvecs(void){
    MatrixXd H = get_H();
    SelfAdjointEigenSolver<MatrixXd> es(H);
    MatrixXd vecs = es.eigenvectors();
    return vecs;
}


 
// test =================================================================


//void test_vector_in_vvector(){}


void test_basis(){
    Hubbard_1d model(4, 2, 1., 6.);
    vector<vector<int>> basis = model.get_basis();
    print4(basis);
    vector<int> bb = {1, 3, 4, 5};
    bool tr = vector_in_vvector(bb, basis);
    cout << tr << endl;
}

void test_double_occ(){
    Hubbard_1d model(6, 3, 1., 6.);
    vector<int> state = {0, 1, 2};
    //cout << model.count_double_occ(state) << endl;
}

void test_hop(){
    /*
    Hubbard_1d model(6, 3, 1., 6.);
    //vector<int> state = {0, 1, 2};
    vector<int> state = {6, 7, 8};
    tuple<vector<int>, int> t = model.hop(state, 2, 1);
    vector<int> state_new = std::get<0>(t);
    int parity = std::get<1>(t);
    print1(state_new);
    cout << parity << endl;
    */
    Hubbard_1d model(2, 1, 1., 6.);
    vector<int> s = {0, 3};
    //print1( model.hop )
 
}


void test_all_hop()
{
    //Hubbard_1d model(6, 3, 1., 6.);
    //vector<int> state = {0, 1, 2, 6, 7, 8};   

    //Hubbard_1d model(4, 2, 1., 6.);
    //vector<int> state = {0, 1, 4, 5};   

    //Hubbard_1d model(8, 4, 1., 6.);
    //vector<int> state = {0, 1, 2, 3, 8, 9, 10, 11};   

    Hubbard_1d model(2, 1, 1., 6.);
    vector<vector<int>> basis = model.get_basis();   

    for (int i = 0; i < basis.size(); i++)
    {
        vector<int> state = basis[i];
        cout << "state" << endl;
        
        print1(state); 
        tuple< vector<vector<int>>, vector<int> > tu;
        tu = model.all_hop(state);
        vector<vector<int>> state_hopped = std::get<0>(tu);
        cout << "state_hopped" << endl;
        print4(state_hopped);
    }
}


void test_get_T(){
    //Hubbard_1d model(6, 3, 1., 6.);
    //vector<int> state = {0, 1, 2, 6, 7, 8};   

    Hubbard_1d model(2, 1, 1., 6.);
    //vector<int> state = {0, 1, 4, 5};   

    MatrixXd T_matr = model.get_T();
    cout << T_matr << endl;
}

void test_get_U(){
    Hubbard_1d model(2, 1, 1., 6.);
    //vector<int> state = {0, 1, 4, 5};   

    MatrixXd U_matr = model.get_U();
    cout << U_matr << endl;
}

void test_get_H(){
    Hubbard_1d model(2, 1, 1., 6.);
    //vector<int> state = {0, 1, 4, 5};   

    MatrixXd H = model.get_H();
    cout << H << endl;
}

void test_eig(){
    Hubbard_1d model(2, 1, 1., 6.);
    //vector<int> state = {0, 1, 4, 5};   

    MatrixXd Es = model.get_eigvals();
    cout << Es << endl;
    MatrixXd vecs = model.get_eigvecs();
    cout << vecs << endl;


}

// main ====================================================

int main(){
    //test_basis();
    //test_hop();
    //test_all_hop();
    //test_get_T();
    //test_get_U();
    test_get_H();
    //test_eig();

}


