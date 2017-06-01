#include <iostream>
#include <iomanip>
#include "data.h"
#include "util.h"

using namespace std;

int main(){
    dataSet test_data;

    parseData("testdata.txt", test_data);
    int m=test_data.get_m();
    int n=test_data.get_n();

    cout<<"m: "<<m<< " n: "<<n<<endl;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            cout<<setw(5)<<test_data.get_dataMtx(i,j);
        }
        cout<<endl;
    }

    

    return 0;
}
