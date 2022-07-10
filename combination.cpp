#include<iostream>
#include<vector>
#include <eigen3/Eigen/Eigen>
using namespace std;
using namespace Eigen;

/*
void all_combines(int a[], int b[], int n, int m, const int M)
{
    if (m>0)
    {
        for (int i=n; i>=m; i--)
        {
            b[m-1] = a[i-1];
            all_combines(a, b, i-1, m-1, M);
        }
    }
    else
    {
        for (int i=0; i < M; i++)
            cout << b[i] << "";
            cout << endl;        
    }
    cout << m << endl;
}
*/

void print1(vector<int> arr)
{
    cout << "begin print arr" << endl;
    for (int i = 0; i < arr.size(); i++) { 
        cout << arr[i] << endl;
    }
    cout << "end print arr" << endl;
}


void print2(vector<vector<int>> arr)
{
    for (int i = 0; i < arr.size(); i++) {
        //cout << arr[i].size() << endl;
        cout << arr[i][0] << endl; 
    }
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
        for (int j=0; j < empty.size(); j++){
            int lastnum = empty[j].back();
            for (int k=(lastnum+1); k <= arr.back(); k++){
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


void vector_to_Vector(vector<int> vec, VectorXd & Vec)
{
    //cout << vec[0] << endl; 
    //Vec[1] = 100;
    for (int i=0; i < vec.size(); i++){
        Vec[i] = vec[i];  
    }
}


int main(){
    /*
    int arr[] = {1,2,3,4};
    int b[4];
    all_combines(arr, b, 4, 3, 3);
    */

    /*    
    vector<vector<int>> empty;
    vector<int> arr = {1, 2, 3, 4};
    comb(empty, arr, 4, 2, 2);
    print2(empty);   
    //cout << empty.size() << endl;   
    */

    /*
    Eigen::VectorXd p(4);
    //cout << p[0] << p[1] << p[2] << p[3] << endl;
    //p << 1, 2, 3, 4;
    //cout << p[0] << p[1] << p[2] << p[3] << endl;
    vector<int> arr = {1, 2, 3, 4};
    vector_to_Vector(arr, p); 
    cout << p[0] << p[1] << p[2] << p[3] << endl;
    */


    
    vector<int> arr = {0, 1, 2, 3};
    vector<vector<int>> brr = get_comb(arr, 2);
    print2(brr);   
    
}



