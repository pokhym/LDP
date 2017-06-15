#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;
int main(){
    int n, r;
    cin>>n;
    cin>>r;

    vector<bool> v(n);
    fill(v.begin(), v.begin()+r, true);

    do{
        for(int i=0; i<n; ++i){
            if(v[i]){
                cout<<(i+1)<<" ";
            }
        }
        cout<<"\n";
    } while(prev_permutation(v.begin(), v.end()));
    return 0;
}
