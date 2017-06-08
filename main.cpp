#include <iostream>
#include <iomanip>
#include <vector>
#include "data.h"
#include "util.h"

using namespace std;

int main(){
    dataSet *test_data=new(dataSet);

    parseData("testdata.txt", *test_data);

    cout<<"parse end"<<endl;
    int m=test_data->get_m();
    int n=test_data->get_n();

    cout<<"m: "<<m<< " n: "<<n<<endl;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            //test_data->set_dataMtx(i+1,j+1,i+j);
            cout<<setw(8)<<" "<<test_data->get_dataMtx(i,j);
        }
        cout<<endl;
    }
    

    delete test_data;
/*
 *    cout<<endl;
 *
 *     double avg_real=0;
 *    int asdf=0;
 *    for(int i=0; i<test_data.get_m(); i++){
 *        if(test_data.get_dataMtx(i,1)!=0){
 *            avg_real=avg_real+test_data.get_dataMtx(i, 1);
 *            asdf++;
 *        }
 *    }
 *    avg_real=avg_real/asdf;
 *
 *    vector<double> outlier={0};
 *    normalizeNeg1toPos1(test_data, outlier);
 *     for(int i=0; i<m; i++){
 *        for(int j=0; j<n; j++){
 *            cout<<setw(8)<<" "<<test_data.get_dataMtx(i,j);
 *        }
 *        cout<<endl;
 *    }
 *
 *
 *    tuplePerturbation(test_data, 3.0);
 *    
 *    for(int i=0; i<m; i++){
 *       for(int j=0; j<n; j++){
 *           cout<<setw(8)<<" "<<test_data.get_dataMtx(i,j);
 *       }
 *       cout<<endl;
 *    }
 *    
 *    double avg=0;
 *    for(int i=0; i<test_data.get_m(); i++){
 *        avg=avg+test_data.get_dataMtx(i, 1);
 *    }
 *    avg=avg/test_data.get_m();
 *
 *
 *    cout<<"avg real: "<<avg_real<<endl;
 *    cout<<"avg est: "<<avg<<endl;
 *
 */
    return 0;
}
