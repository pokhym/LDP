#include "util.h"
#include <random>
#include <iomanip> // remove later
#include <cstdlib>
#include <ctime>

#define TAB 0x09
#define LARGENUMBER 999999999

using namespace std;

/* parseData
 * DESCRIPTION: Parses a comma separated list of data. Each column represents
 * an attribute an each row represents a user
 * INPUTS: filename, dataset
 * OUTPUTS: 0 on success -1 on failure
 */
int parseData(char const* filename, dataSet &data, int init_m, int init_n){
    if(init_m<1 || init_n<1 || filename==NULL)
        return -1;

    // cout<<"initial resize"<<endl;
    data.rekeep_dataMtx(init_m,init_n);
    // cout<<"initial resize complete"<<endl;

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
            // if(row_count%10000==0)
                // cout<<"row count: "<<row_count<<endl;
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
    else
        return -1;

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
                    data.set_dataMtx(m+1,n+1, LARGENUMBER);
                }
            }
        }
    }

    pair<vector<double>, vector<double>> max_min_pair;
    max_min_pair.first=max;
    max_min_pair.second=min;

    return max_min_pair;
}

/* factorial
 * DESCRIPTION: Calcualtes the factorial of the input with the tgamma function
 * INPUTS: long long representing the integer factorial we want to compute
 * OUTPUTS: the long long result of i!
 */
long long factorial(long long i){
    if(i==1)
        return 1;
    else
        return i*factorial(i-1);
}

/* nCr
 * DESCRIPTION: Calculates the nCr combinations of two inputs
 * INPUTS: n, r as integers
 * OUTPUTS: number of combinations
 */
unsigned long long nCr(long long n, long long r){
    // use the gamma function to efficiently calculate the factorial funtion
    cout<<"nCr func: "<< tgamma(n+1)/(tgamma(r+1)*tgamma(n-r+1))<<endl;

    return tgamma(n+1)/(tgamma(r+1)*tgamma(n-r+1));
}

/* tuplePerturbation
 * DESCRIPTION: Implements algorithm 2 in the harmony algorithm
 * INPUTS: dataSet containing all the columns, epsilon, and vector of outlier
 * conditions
 * OUTPUTS: perturbed value on success -1 on fail
 */
vector<double> tuplePerturbation(dataSet &data, double epsilon){
    if(data.get_n()==0 || data.get_m()==0 || epsilon<=0)
        return vector<double>(-1);

    vector<double> perturbed(data.get_m(), 0.0); // store result here
    int d=data.get_n(); // number of attributes

    // initialize rand
    srand(time(0));

    // initialize the random engine for the rest of the distributions
    default_random_engine gen;

    //    --
    //    | 2^(d-1), if d is odd
    // Cd=|
    //    | 2^(d-1)- 1/2 * nCr(d, d/2), otherwise
    //    --
    //   --
    //   |      2^d + Cd*(exp(epsilon) - 1)
    //   | ------------------------------------ , if d is odd
    //   | nCr(d-1, (d-1)/2)*(exp(epsilon) - 1)
    // B=|
    //   |    2^d + Cd*(exp(epsilon) - 1)
    //   | -------------------------------- , otherwise
    //   | nCr(d-1, d/2)*(exp(epsilon) - 1)
    //   --
    // First calculate B
    long double B=0.0;
    if(d%2==0){

        // 2^d     + Cd                *                                             exp(epsilon)-1
        B=(pow(2, d)+(pow(2, d-1)-0.5*(tgamma(d+1)/(tgamma(d/2+1)*tgamma(d-d/2+1))))*(exp(epsilon)-1))/
          ((tgamma(d-1+1)/(tgamma((d-1)/2+1)*tgamma(d-1-(d-1)/2+1)))*(exp(epsilon)-1));
           // ncR(d-1, d/2)                   *                      exp(epsilon)-1
    }
    else{
        B=(pow(2, d)+pow(2, d-1)*(exp(epsilon)-1))/
          ((tgamma(d-1+1)/(tgamma(d/2+1)*tgamma(d-1-d/2+1)))*(exp(epsilon)-1));
    }

    // for each row
    for(int i=0; i<data.get_m(); i++){
        /////////////////////////////////// STEP 1
        vector<int> v(data.get_n(), 0);
        // for each column generate v
        for(int j=0; j<data.get_n(); j++){
            // grab Aj a random attribute
            // double Aj=rand()%d; // needs to have a max number of d-1 attributes (0 indexed)

            // generate a v vector containing {1,-1} for each Aj with following prob dist
            // NOTE: Aj=j in this case
            //             --
            //             | 1/2 + 1/2 * ti[Aj], if x=1
            // Pr[v[Aj]=x]=|
            //             | 1/2 - 1/2 * ti[Aj], if x=-1
            //             --
            // 0.5+0.5*1=.75 0.5+0.5*.75=.875
            // 0.5-0.5*1=.25 0.5-0.5*.75=.125
            //-------------------------------
            //           1.0             1.0 always sums to 1 so only need to use one of the lines
            // only do this if not an outlier
            if(data.get_dataMtx(i,j)!=LARGENUMBER){
                double prob_v=0.5+0.5*data.get_dataMtx(i, j);
                bernoulli_distribution distty(prob_v); // use the result of 1 to mark as 1 else -1

                if(distty(gen)){
                    v[j]=1;
                }
                else{
                    v[j]=-1;
                }
            }
            else{ // otherwise set LARGENUMBER use as a flag for later
                v[j]=LARGENUMBER;
            }
        }

        //////////////////////////////// STEP 2
        // Time to choose whether we want to use T+ or T-
        bernoulli_distribution dist(exp(epsilon)/(exp(epsilon)+1));
        int T_plus_minus=dist(gen);

        // Now generate T+ and T-
        // T+ is defined as t* element {B-, B+}^d where t* * v>0
        // T- is defined as t* element {B-, B+}^d where t* * v<=0
        // For example from the paper if ti=< 1, 1 >
        // T+ = { <B, B> }
        // T- = { <-B, -B>, <-B, B>, <B, -B> }
        // Instead of choosing one combination from many we just generate one at random
        // For example in the case of T+ the idea is to generate the indexes we want to be positive
        // then fill in the gaps with the rest as negative

        // Case T+: There must be more than d/2 positive B's
        if(T_plus_minus){
            int num_positive=rand()%(d/2)+d/2+1; // if d=20 we want to have between 11-20 B
                                                 // achieve by doing d/2+1=11 and rand()%(d/2)=0-9

            vector<double> B_flags(d, 0.0);
            // Set the positive values first
            // there are 2 cases accounting for v=-1 and v=1
            int count=0;
            while(count!=num_positive){
                int idx=rand()%d;
                if(v[idx]==1 && B_flags[idx]==0.0){ // set 1 for positive vars
                    B_flags[idx]=1;
                    count++;
                }
                else if(v[idx]==-1 && B_flags[idx]==0.0){ // v=-1 and we need to set -1
                    B_flags[idx]=-1;
                    count++;
                }
                else{
                    count++;
                }
                // ignore LARGENUMBER case
            }

            // we can now set the negative values
            for(int p=0; p<(int)B_flags.size(); p++){
                if(B_flags[p]==0){
                    if(v[p]==1){
                        B_flags[p]=-1;
                    }
                    else if(v[p]==-1){
                        B_flags[p]=1;
                    }
                    // ignore LARGENUMBER case
                }
            }

            // we have now obtained the correct signs of B in the current user row
            // set the values accordingly
            double sum=0;

            for(int p=0; p<d; p++){
                if(B_flags[p]!=0.0){
                    data.set_dataMtx(i+1, p+1, B*B_flags[p]);
                    sum=sum+(double)v[p]*data.get_dataMtx(i, p);
                }
            }

            // if(sum<=0){
            //     cout<<"T+ i "<<i<<" sum "<<sum<<" "<<B_flags.size()<<endl;
            //     for(int i=0; i<(int)B_flags.size(); i++){
            //         cout<<B_flags[i]<<" "<<v[i]<<" "<<" ***** "<<B_flags[i]*(double)v[i]<<" **** ";
            //     }
            //     cout<<endl; cout<<endl;
            // }
        }

        // Case T-: There must be less than or equal to d/2 negative B's
        else{
            int num_positive=rand()%(d/2+1); // if d=20 we want to have between 0-10 B
                                             // achieve by doing rand()%(d/2+1)=0-10
            vector<double> B_flags(d, 0.0);
            // Set the positive values first
            // there are 2 cases accounting for v=-1 and v=1
            int count=0;
            while(count!=num_positive){
                int idx=rand()%d;
                if(v[idx]==1 && B_flags[idx]==0.0){ // set 1 for positive vars
                    B_flags[idx]=1;
                    count++;
                }
                else if(v[idx]==-1 && B_flags[idx]==0.0){ // v=-1 and we need to set -1
                    B_flags[idx]=-1;
                    count++;
                }
                else{
                    count++;
                }
                // ignore LARGENUMBER case
            }

            // we can now set the negative values
            for(int p=0; p<(int)B_flags.size(); p++){
                if(B_flags[p]==0){
                    if(v[p]==1){
                        B_flags[p]=-1;
                    }
                    else if(v[p]==-1){
                        B_flags[p]=1;
                    }
                }
            }

            // we have now obtained the correct signs of B in the current user row
            // set the values accordingly
            double sum=0;
            for(int p=0; p<d; p++){
                if(B_flags[p]!=0){
                    data.set_dataMtx(i+1, p+1, B*B_flags[p]);
                    sum=sum+(double)v[p]*data.get_dataMtx(i, p);
                }
            }
            // if(sum>0){
            //     cout<<"T- i "<<i<<" sum "<<sum<<endl;
            // }
        }

    } // end each row

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
