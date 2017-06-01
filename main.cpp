#include <iostream>
#include "data.h"
#include "util.h"

using namespace std;

int main(){
    dataSet test_data;

    test_data.rekeep_dataMtx(3,2);
    for(int i=0; i<2; i++){
        for(int j=0; j<3; j++){
            test_data.set_dataMtx(i,j,i+j);
        }
    }

    vector<double> asdf=tuplePerturbation(test_data, 0.5);
    
/*
 *    parseData("testdata.txt", test_data);
 *    
 *    int m=test_data.get_m();
 *    int n=test_data.get_n();
 *
 *    cout<<"m: "<<m<< " n: "<<n<<endl;
 *
 *    for(int i=0; i<n; i++){
 *        for(int j=0; j<m; j++){
 *            cout<<test_data.get_dataMtx(i,j)<<endl;
 *        }
 *    }
 */
    return 0;
}
