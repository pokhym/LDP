#include "util.h"
#include <random>
#include <iomanip> // remove later

#define TAB 0x09
#define LARGENUMBER 9999999999

using namespace std;

/* parseData
 * DESCRIPTION: Parses a comma separated list of data. Each column represents
 * an attribute an each row represents a user
 * INPUTS: filename, dataset
 * OUTPUTS: 0 on success -1 on failure
 */
int parseData(char const* filename, dataSet &data){
    cout<<"initial resize"<<endl;
    data.rekeep_dataMtx(1000000,100);
    cout<<"initial resize complete"<<endl;

    string line;
    ifstream myfile(filename);
    int prev_data_idx=0; // mark the end of the previous data
    int row_count=0;
    int col_count=0;
    int max_col=0;
    double val=0;

    // open file
    if(myfile.is_open()){
        // read all lines
        while(getline(myfile, line)){
            row_count++;
            if(row_count%10000==0)
                cout<<"row count: "<<row_count<<endl;
            for(int i=0; i<(int)line.length(); i++){
                // check for the delimiter
                // parse the number from the last idx
                // 2 cases
                // 1. ,0.3,
                // 2. ,0.3\n
                // NOTE: prev_data_idx and i are both representative of the end of
                // the string aka it is either , or \n
                if(line[i]==','){
                    col_count++;
                    if(col_count>max_col){
                        max_col=col_count;
                    }

                    // weird case
                    if(prev_data_idx!=0){prev_data_idx++;}

                    // enter null terminated
                    char temp_num_char[i-prev_data_idx+1];
                    temp_num_char[i-prev_data_idx]='\0';

                    // copy
                    for(int j=0; j<i-prev_data_idx; j++){
                        temp_num_char[j]=line[j+prev_data_idx];
                    }

                    // update next string beginning position
                    prev_data_idx=i;

                    // convert into a string so we can use stod
                    string temp_num_string=temp_num_char;
                    val=stod(temp_num_string);

                    data.rekeep_dataMtx(row_count, max_col);
                    data.set_dataMtx(row_count, col_count, val);
                }

                else if(line[i+1]=='\0'){
                    col_count++;
                    if(col_count>max_col){
                        max_col=col_count;
                    }

                    if(prev_data_idx!=0){prev_data_idx++;}

                    char temp_num_char[i-prev_data_idx];
                    temp_num_char[i-prev_data_idx+1]='\0';
                    for(int j=0; j<=i-prev_data_idx; j++){
                        temp_num_char[j]=line[j+prev_data_idx];
                    }
                    prev_data_idx=i;

                    string temp_num_string=temp_num_char;
                    val=stod(temp_num_string);

                    data.rekeep_dataMtx(row_count, max_col);
                    data.set_dataMtx(row_count, col_count, val);
                }
            }

            prev_data_idx=0;
            col_count=0;
        }
    }

    return 0;
}

/* normalizeNeg1toPos1
 * DESCRIPTION: normalizes each column of the input data set -1 to 1
 * INPUTS: dataSet, outlier condition
 * OUTPUTS: 0 on success -1 on fail
 */
pair<vector<double>, vector<double>> normalizeNeg1toPos1(dataSet &data, vector<double> outlier){
    //if(data.get_m()==1 || data.get_n()==1 || outlier.size()==0)
        //return pair<vector<double>, vector<double>>(0);

    // store the max values of each for normalization
    vector<double> max(data.get_n(), -LARGENUMBER);
    vector<double> min(data.get_n(), LARGENUMBER);
    vector<double> scales(data.get_n(), 0.0);

    // loop through rows and columns to obtain max
    double max_curr_col=-LARGENUMBER; // maximum in current column
    double min_curr_col=LARGENUMBER; // minimum in current column
    for(int n=0; n<data.get_n(); n++){
        for(int m=0; m<data.get_m(); m++){
            if(data.get_dataMtx(m,n)!=outlier[0]){ // hard code to 0 outlier
                if(data.get_dataMtx(m,n)>max_curr_col)
                    max_curr_col=data.get_dataMtx(m,n); // update

                if(data.get_dataMtx(m,n)<min_curr_col)
                    min_curr_col=data.get_dataMtx(m,n); // update
            }
        }
        // save into vector and iterate to next column
        max[n]=max_curr_col;
        min[n]=min_curr_col;

        // reset min_cur_col
        max_curr_col=-LARGENUMBER;
        min_curr_col=LARGENUMBER;
    }


    // normalize per column
    for(int n=0; n<data.get_n(); n++){
        for(int m=0; m<data.get_m(); m++){
            // if maximum is 0 (all outliers) don't do anything
            if(max[n]==0){}

            // else we need to normalize all nonzero values
            else{
                // if the data in that m,n is not an outlier normalize
                if(data.get_dataMtx(m,n)!=outlier[0]){ // hardcode outlier 0
                    double new_val=data.get_dataMtx(m,n);
                    new_val=(new_val-min[n])/(max[n]-min[n]);

                    // center around -1 to 1
                    new_val=new_val-0.5;
                    new_val=new_val*2.0;

                    // save scales
                    scales[n]=scales[n]+new_val;

                    // save new value into data
                    data.set_dataMtx(m+1,n+1,new_val);
                }

                // else we need to set to 0
                else{
                    data.set_dataMtx(m+1,n+1,0.0);
                }
            }
        }
    }

    pair<vector<double>, vector<double>> max_min_pair;
    max_min_pair.first=max;
    max_min_pair.second=min;

    return max_min_pair;
}


/* tuplePerturbation
 * DESCRIPTION: Implements algorithm 2 in the harmony algorithm
 * INPUTS: dataSet containing all the columns, epsilon, and vector of outlier
 * conditions
 * OUTPUTS: perturbed value on success -1 on fail
 */
vector<double> tuplePerturbation(dataSet &data, double epsilon){
    if(data.get_n()==0 || data.get_m()==0)
        return vector<double>(-1);

    vector<double> perturbed(data.get_m(), 0.0); // store result here

    // uniformally decide which attribute to calculate
    // currently only one so ignore this step
    default_random_engine generator;
    uniform_int_distribution<int> distribution(0,data.get_n()-1);
    random_device rdd;
    mt19937 genn(rdd());



    // loop through all the different entries (not attributes)
    for(int i=0; i<data.get_m(); i++){
        double probability, numerator, denominator;
        int n=distribution(generator);
        //          ti[Aj]*(exp^(epsilon)-1)+exp^(epsilon)+1
        // Pr[u=1]= ----------------------------------------
        //                     2*exp^(epsilon)+2

        if(data.get_dataMtx(i,n)!=0){
            numerator=data.get_dataMtx(i, n)*(exp(epsilon)-1.0)+exp(epsilon)+1.0;
            denominator=(double)2.0*exp(epsilon)+2.0;
            probability=numerator/denominator;

            // generate random crap code pulled crom c++reference
            bernoulli_distribution d(probability); // set p

            // grab data values from bernoulli
            int bernoulli_result=d(genn);

            // assign values
            if(bernoulli_result){
                perturbed[i]=(double)(exp(epsilon)+1.0)/(exp(epsilon)-1.0);//(data.get_n())*(exp(epsilon)+1.0)/(exp(epsilon)-1.0);
                data.set_dataMtx(i+1, n+1, (double)(exp(epsilon)+1.0)/(exp(epsilon)-1.0));//(data.get_n())*(exp(epsilon)+1.0)/(exp(epsilon)-1.0));
            }
            else{
                perturbed[i]=(double)-1.0*(exp(epsilon)+1.0)/(exp(epsilon)-1.0);//-1.0*(data.get_n())*(exp(epsilon)+1.0)/(exp(epsilon)-1.0);
                data.set_dataMtx(i+1, n+1, (double)-1.0*(exp(epsilon)+1.0)/(exp(epsilon)-1.0));//(double)-1.0*(data.get_n())*(exp(epsilon)+1.0)/(exp(epsilon)-1.0));
            }

            for(int idx=0; idx<data.get_n(); idx++){
                if(idx==n){}
                else{
                    data.set_dataMtx(i+1, idx+1, 0.0);
                }
            }
        }
    }

    return perturbed;
}

/* undoNorm
 * DESCRIPTION: Undos the [-1, 1] normalization
 * INPUTS: max_min_pair with first=max, second=min, and a preNormalization
 * vector of values
 * OUTPUTS: 0 on success -1 on failure
 */
int undoNorm(std::pair<std::vector<double>, std::vector<double>> max_min_pair, std::vector<double> &preNorm){
    for(int i=0; i<(int)preNorm.size(); i++){
        preNorm[i]=((preNorm[i]+1)/2)*(max_min_pair.first[i]-max_min_pair.second[i])+max_min_pair.second[i];
    }
    return 0;
}
