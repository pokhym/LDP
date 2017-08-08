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

/* gen_rand_mtx
 * DESCRIPTION: used in bassily_smith to generate a random binary mtx
 * INPUTS: dimensions m x k, m=number of bits per message, k=number of unique msgs
 * OUTPUTS: none, this function will edit a preallocated dataset
 */
 void gen_rand_mtx(boolSet &rand_mtx, int m, int k){
     static mt19937_64 gen(time(0));
     static uniform_int_distribution<int> _0_1(0,1);
     int total=0; int plus=0;
     double sum=0;

     // randomly generate plus or minus 1/sqrt(m) for each
     cout<<rand_mtx.get_m()<<" "<<rand_mtx.get_n()<<endl;
     long long max=rand_mtx.get_m()*rand_mtx.get_n();
    // long long i=0;

     // weird error case with 0,0 generate artificially first this may be a
     // liability for security but w/e
    //  if(_0_1(gen)){rand_mtx.set_dataMtx(0+1, 0+1, 1.0/sqrt(m));}
    //  else{rand_mtx.set_dataMtx(0+1, 0+1, -1.0/sqrt(m));}
     //
    //  while(i<(max-1)){
    //      int _m=rand()%m;
    //      int _k=rand()%k;
    //      if(_0_1(gen)){
    //          if(rand_mtx.get_dataMtx(_m, _k)==0){
    //              rand_mtx.set_dataMtx(_m+1, _k+1, (double)1.0/sqrt(m));
    //              total++;
    //              plus++;
    //              sum+=(double)1.0/sqrt(m);
    //              i++;
    //          }
    //      }
    //      // 0 case: -1/sqrt(m)
    //      else{
    //          if(rand_mtx.get_dataMtx(_m, _k)==0){
    //              rand_mtx.set_dataMtx(_m+1, _k+1, (double)-1.0/sqrt(m));
    //              total++;
    //              sum+=-(double)1.0/sqrt(m);
    //              i++;
    //          }
    //      }
    //      if(i%10000==0){cout<<"total: "<<max<<" "<<(double)100*i/max<<" % \r"; cout.flush();}
    //  }
    //  sum=0;

    int j=0;
    int i=0;
    double sumsum=0;
    for(i=0; i<m; i++){
        for(j=0; j<k; j++){
            if(_0_1(gen)){
                rand_mtx.set_boolMtx(i+1,j+1,1);
            }
            else{
                rand_mtx.set_boolMtx(i+1,j+1,0);
            }
        }
    }

    // calculate inner product
    // long long j=0;
    // double sumsum=0;
    // int asdf=0;
    // for(asdf=0; asdf<k; asdf++){
    //     for(j=0; j<k; j++){
    //         for(i=0; i<m; i++){
    //             if(asdf!=j){
    //                 sum+=rand_mtx.get_dataMtx(i, asdf)*rand_mtx.get_dataMtx(i, j);
    //             }
    //         }
    //         //cout<<j<<" "<<sum<<endl;
    //         sumsum+=sum;
    //         cout<<"non self inner product: "<<sumsum/k<<endl;
    //         sum=0;
    //     }
    //     sumsum=0;
    // }
    // sumsum=0; sum=0;
    //
    // for(j=0; j<k; j++){
    //     for(i=0; i<m; i++){
    //         sum+=rand_mtx.get_dataMtx(i, j)*rand_mtx.get_dataMtx(i, j);
    //     }
    //     //cout<<j<<" "<<sum<<endl;
    //     sumsum+=sum;
    //     sum=0;
    // }
    cout<<"self inner product: "<<sumsum/k<<endl;
     //cout<<"percent + generated: "<<(double)plus/total<<" sum: "<<sum<<endl;

     // verify inner product of columns with itself seems to be working
     // but then this doesn't account for the fact that he difference maybe less than
     // machine epsilon?  it says sum!=1 when it prints sum=1

 }

/* bassily_smith
 * DESCRIPTION: bassily_smith implementation of frequency estimate
 * INPUTS: a data set, epsilon (privacy budget), beta (confidence in error bound)
 * OUTPUTS: k number of frequncy estimates (k is the number of unique data vals)
 */
 vector<double> bassily_smith(dataSet &data, double epsilon, double beta, int k){
     double gamma, m;

     // calculate gamma (some crap used to determine number of bits per msg)
     //                log2(2k/beta)
     // gamma = sqrt(-----------------)
     //                 epsilon^2*n
     // numerator=log2(2*k/beta);
     // denominator=epsilon*epsilon;
     gamma=(double)sqrt((log(2.0*k/beta))/(epsilon*epsilon*data.get_m())); //sqrt(numerator/denominator);

     // calculate m (number of bits per msg)
     //      log2(k+1)*log2(2/beta)
     // m = ------------------------
     //            gamma^2
     // numerator=log2(k+1)*log2(2/beta);
     // denominator=gamma*gamma;
     m=(double)(log(k+1.0)*log(2.0/beta))/(gamma*gamma); //numerator/denominator;
     cout<<m<<endl;
     if((m-(int)m)!=0){m=(int)m+1;} // update extra m, eg m=4.1 -> m=5 needs to be whole

     // generate random matrix
     boolSet phi;
     cout<<"m: "<<m<<" k: "<<k<<endl;
     phi.rekeep_boolMtx(m, k);
     //(phi, m, k);

     // we will now create the z_i vectors (size k).  for each user
     // these contain one data value and the rest are zerO
     mt19937_64 gen(time(0));
     uniform_int_distribution<int> idx_chooser(0,m-1);

     mt19937_64 gen1(time(0));
     bernoulli_distribution _0_1((double)exp(epsilon)/(exp(epsilon)+1.0));

     vector<pair<int, double> > zi_holder;

     cout<<"for each user"<<endl;

     for(int i=0; i<data.get_m(); i++){
         // draw a uniform index
         int idx=idx_chooser(gen);

         // draw a 0 or 1 based on the following probability
         //             exp(expsilon)
         // prob = ---------------------
         //          exp(expsilon) + 1
         //int res=_0_1(gen1);
         //if(res==0){cout<<res<<endl;}

         // let c_epsilon equal the below
         //              exp(expsilon+1)
         // c_epsilon = ------------------
         //              exp(expsilon-1)
         double c_epsilon=(double)(exp(epsilon)+1.0)/(exp(epsilon)-1.0);

         // if 1 assign value alpha to z_i where alpha is equal to the following
         // alpha = c_epsilon * m * phi[s, ti[Aj]]
         double alpha;
         if(_0_1(gen1)){ // NOTE: assume data starts indexing at 1
             double actual_value;
             if(phi.get_boolMtx(idx, data.get_dataMtx(i,0)-1)==1){actual_value=1.0/sqrt(m);}
             else{actual_value=-1.0/sqrt(m);}
             alpha=(double)c_epsilon*m*actual_value;
         }

         // else 0 assign value alpha to z_i where alpha is equal to the following
         // alpha = -c_epsilon * m * phi[s, ti[Aj]]
         else{
             double actual_value;
             if(phi.get_boolMtx(idx, data.get_dataMtx(i,0)-1)==1){actual_value=1.0/sqrt(m);}
             else{actual_value=-1.0/sqrt(m);}
             alpha=(double)-c_epsilon*m*actual_value;
         }

         // store the pair value per user as zi=pair<int, double>
         zi_holder.push_back(pair<int, double>(idx, alpha));
     }

     cout<<"calculate zbar"<<endl;

     // calculate z_bar average over n users (data.get_m)
     vector<double> z_bar(m, 0.0); cout<<zi_holder.size()<<" "<<data.get_m()<<endl;

     for(int i=0; i<(int)data.get_m(); i++){
         // cout<<zi_holder[i].first<<" ";
         //cout<<zi_holder[i].first<<endl;
         z_bar[zi_holder[i].first]=z_bar[zi_holder[i].first]+zi_holder[i].second;
         //cout<<zi_holder[i].first<<" "<<zi_holder[i].second<<endl;
     }
     //cout<<endl;
     for(int i=0; i<m; i++){
         z_bar[i]=(double)z_bar[i]/data.get_m();
         //cout<<z_bar[i]<<" ";
     }
     //cout<<endl;

     cout<<"calculate k_freq_ests"<<endl;

     // for each unique value estimate the frequency by obatining the innerproduct
     // between each column of phi and z_bar
     vector<double> k_freq_ests(k, 0.0);
     for(int j=0; j<k; j++){
         for(int i=0; i<phi.get_m(); i++){
             k_freq_ests[j]=k_freq_ests[j]+z_bar[i]*phi.get_boolMtx(i, j)*;
         }
         //cout<<k_freq_ests[j]<<" ";
     }
     //cout<<endl;

     return k_freq_ests;
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
