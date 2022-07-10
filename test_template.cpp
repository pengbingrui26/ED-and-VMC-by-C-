# include<iostream>
# include<vector>
using namespace std;

template<typename T>
void myprint(T vec)
{
    for (int i=0; i < vec.size(); i++)
    {
        cout << vec[i] << " ";
    }
    cout << endl;
}


template<typename T1, typename T2>
void myprint_two(T1 vec1, T2 vec2)
{
    for (int i=0; i < vec1.size(); i++)
    {
        cout << vec1[i] << " ";
    }
    cout << endl;
    for (int i=0; i < vec2.size(); i++)
    {
        cout << vec2[i] << " ";
    }
    cout << endl;
}


int main()
{
    //vector<int> v = {1,2,3};
    //myprint<vector<int>>(v);

    vector<int> v1 = {1,2,3};
    vector<double> v2 = {11.,12.,13., 14.};
    myprint_two<vector<int>, vector<double>>(v1, v2);

}




