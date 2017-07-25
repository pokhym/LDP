#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "data.h"
#include "util.h"

using namespace std;
#define LARGENUMBER 999999999

void printall(dataSet &data){
    // int m=data.get_m();
    // int n=data.get_n();
    //
    // for(int i=0; i<m; i++){
    //     for(int j=0; j<n; j++){
    //         cout<<setw(12)<<data.get_dataMtx(i,j);
    //     }
    //     cout<<endl;
    //     cout<<"next row"<<endl;
    // }
}

double MSE(vector<double> avg_real, vector<double> avg_est){
    int n=avg_real.size();
    double total=0.0;
    for(int i=0; i<n; i++){
        cout<<"total MSE: "<<total<<endl;
        total=total+(avg_real[i]-avg_est[i])*(avg_real[i]-avg_est[i]);
    }
    return (double)total/n;
}


int main(int argc, const char* argv[]){
    if(argc<4 || argc>5){
        cout<<"usage: ./test path_to_data initial_m initial_n"<<endl;
	cout<<"usage: ./test path_to_data initial_m initial_n epsilon"<<endl;
        cout<<argc<<endl;
        return 0;
    }
    double sum=0.0;
    double avg_mse=0.0;
    int m;
    int n;
    int rep;
    double epsilon=0.1;

    // case for one round
    if(argc==5){
        rep=0;
    }
    else
	rep=7;

    for(int o=0; o<(rep); o++){
	if(o==0){epsilon=0.05;}
	else if(o==1){epsilon=0.1;}
	else if(o==2){epsilon=0.2;}
	else if(o==3){epsilon=0.4;}
	else if(o==4){epsilon=0.8;}
	else if(o==5){epsilon=1.6;}
	else if(o==6){epsilon=3.2;}

        for(int p=0; p<20; p++){
            dataSet *test_data=new(dataSet);

            // cout<<"begin parsing data"<<endl;

            // parse data
            int lol;
            // cout<<stoi(argv[2])<<endl;
            // cout<<stoi(argv[3])<<endl;
            lol=parseData(argv[1], *test_data, stoi(argv[2]), stoi(argv[3]));
            if(lol==-1){
                cout<<"failed to parse data you might want to check the path"<<endl;
                return -1;
            }

            printall(*test_data);

            // cout<<"parse end"<<endl;

            m=test_data->get_m();
            n=test_data->get_n();
            // cout<<"m (users): "<<m<<" n (attributes): "<<n<<endl;

            // calculate real averages
            vector<double> avg_real(test_data->get_n(), 0.0);
            vector<int> asdf(test_data->get_n(), 0);

            for(int i=0; i<test_data->get_m(); i++){
                for(int j=0; j<test_data->get_n(); j++){
                    if(test_data->get_dataMtx(i,j)!=0){
                        avg_real[j]=avg_real[j]+test_data->get_dataMtx(i, j);
                        asdf[j]++;
                    }
                }
            }

            for(int i=0; i<test_data->get_n(); i++){
                avg_real[i]=avg_real[i]/asdf[i];
            }

            // data normalization
            // cout<<"normalizing columns"<<endl;
            vector<double> outlier={0};
            pair<vector<double>, vector<double>> scales=normalizeNeg1toPos1(*test_data, outlier);

            // cout<<"done normalizing"<<endl;

            printall(*test_data);
            // cout<<endl;

            // perturb data
            // cout<<"perturb"<<endl;
            vector<double> perturb_scales=tuplePerturbation(*test_data,epsilon);
            printall(*test_data);
            // cout<<"end perturb"<<endl;

            // calculate estimated average
            vector<double> avg(test_data->get_n(), 0.0);
            vector<int> counters(test_data->get_n(), 0);

            for(int j=0; j<test_data->get_n(); j++){
                for(int i=0; i<test_data->get_m(); i++){
                    if(test_data->get_dataMtx(i,j)!=LARGENUMBER){
                        avg[j]=avg[j]+test_data->get_dataMtx(i,j);
                        counters[j]++;
                    }
                }
            }

            for(int i=0; i<test_data->get_n(); i++){
                avg[i]=avg[i]/counters[i];
            }

            // undo normalizaion
            undoNorm(scales, avg);

            // print real averages and estimated averages
            double rel_error=0.0;
            for(int i=0; i<test_data->get_n(); i++){
                 // cout<<"avg_real: "<<avg_real[i]<<setw(10)<<" avg_est: "<<avg[i]<<endl;
                 rel_error+=abs((avg[i]-avg_real[i])/avg_real[i]);
            }
            // cout<<"percent releative error: "<<rel_error/test_data->get_n()<<endl;
            // cout<<endl;

            double mse=MSE(avg_real, avg);
            //cout<<"MSE: "<<mse<<endl;

            sum=sum+mse;
            avg_mse=sum/(p+1);
            //cout<<"running average("<<p+1<<"): "<<avg_mse<<endl;

            delete test_data;
        }
        cout<<endl;
        cout<<"file name: "<<argv[1]<<" epsilon: "<<epsilon<<endl;
        cout<<"m: "<<m<<" n: "<<n<<" final MSE: "<<avg_mse<<endl;
        cout<<endl;

        sum=0.0; avg_mse=0.0;

        epsilon+=0.1;

    }
    return 0;
}
