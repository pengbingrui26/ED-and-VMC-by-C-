# include<iostream>
# include<vector>
# include<algorithm>
using namespace std;


void print(vector<int> arr)
{
    for (int i = 0; i < arr.size(); i++) {
        cout << arr[i] << endl;
    }
}

void print2(vector<int> arr)
{
    for (int i = 0; i < arr.size(); i++) {
        cout << arr[i] << endl;
    }
}


void increase(vector<int>& arr)
{
    for (int i = 0; i < arr.size(); i++) {
        //cout << i << endl;
        //cout << arr[i] << endl;
        arr[i] ++;        
    }
}


void add(vector<vector<int> >& arr, vector<int>& brr)
{
    for (int i=0; i < brr.size(); i++){
        vector<int> tmp = {brr[i]};
        arr.push_back(tmp);
    }    
}


void new_arr(vector<int>& arr)
{
    vector<int> tmp = {1, 2, 3};
    arr = tmp;
}

vector<int> fun(int N){
    vector<int> vec = {1,2,3};
    vec[0] = N;
    return vec;
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




int main()
{
    //vector<int> a = {1,2,3,4};
    //print(a);
    //increase(a);
    //vector<vector<int>> empty;
    //add(empty, a);
    //cout << "\n" << endl;
    //print(empty);
    //cout << empty.back()[0] << endl;
    //cout << empty.size() << endl;
    //new_arr(a);
    //print(a); 
    //int N = 6;
    //vector<int>v = fun(6);
    //cout << v[0] << endl;


    /*
    vector<int> a = {1,2,3};
    vector<int> b = {1,2,3};
    vector<int> c = combine_two(a, b);
    //print(c); 
    bool tr = (a==b); 
    //cout << tr << endl;   
    for (int i=0; i < 3; i++){
        a[i] -= 5;
    }
    //cout << a[0] << a[1] << a[2] << endl;
    vector<int> p;
    cout << p.size() << endl;
    vector<vector<int>> v;
    v.push_back(p);
    v.push_back(a);
    cout << v.size() << endl;
    */


    /*
    vector<int> a = {1,100,2,3};
    int key = 3;
    
    vector<int>::iterator it = find( a.begin(), a.end(), key);
    //cout << *it << endl;
    cout << std::distance(a.begin(), it) << endl;
    */

    vector<vector<int>> vv;

    vector<int> a = {1,2,3};
    vector<int> b = {4,5,6};
    vector<int> c = {7,8,9};
 
    vv.push_back(a);
    vv.push_back(b);
    vv.push_back(c);
 
    vector<int> key = {7,8,9};  
    vector<vector<int>>::iterator it = find( vv.begin(), vv.end(), key);
    //cout << *it << endl;
    cout << std::distance(vv.begin(), it) << endl;

}


