#include<iostream>
#include<vector>
using namespace std;

void print1(vector<int> arr){
    for (int k = 0; k < arr.size(); k++)
    { 
        cout << arr[k] << endl;
    }
}


void print2(vector<double> arr){
    for (int k = 0; k < arr.size(); k++)
    { 
        cout << arr[k] << endl;
    }
}

void bubble_sort(vector<double>& arr){
    const int N = arr.size(); 
    for (int i=N-1; i>=0; i--)
    {
        for (int j=0; j<=i-1; j++)
        {
            if (arr[j] > arr[j+1]){
                double tmp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp;  
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
                double tmp1 = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp1;

                double tmp2 = idx[j];
                idx[j] = idx[j+1];
                idx[j+1] = tmp2;
            }
        }
    }   
    return idx;
}

int main(){
    vector<double> v = {30, 10., 4.5, 2};
    //bubble_sort(v);
    print2(v);
    print1(argsort(v)); 

}

