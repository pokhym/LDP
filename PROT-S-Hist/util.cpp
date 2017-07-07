#include "util.h"
#include <random>
#include <iomanip> // remove later
#include <cstdlib>
#include <ctime>
#include <bitset>

#define TAB 0x09
#define LARGENUMBER 9999999999
int m_code=31; //  m is the number of bits the randomizer will produce

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

    cout<<"initial resize"<<endl;
    data.rekeep_dataMtx(init_m,init_n);
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

/* code
 * DESCRIPTION: Encodes a thing into binary hamming code
 * INPUT: A number to encode into binary hamming code
 * OUTPUT: A character array of the encoded hamming stuff
 */
 vector<int> code(int input){
     // create binary input, figure out the length of it
     // create hamming code as required

     bitset<32> bit_rep(input); // displays input as binary representation
     int right, left; // 01100 read right to left so init, right=1 left=0
                      // then right=1 left=1 ---->
     int idx_info=31;
     while(1){
         right=bit_rep[idx_info-1];
         left=bit_rep[idx_info];

         if(right==1 && left==0){
             cout<<"input: "<<input<<endl;
             idx_info=idx_info-1; // this is now the number of bits
             cout<<"number of digits: "<<idx_info+1<<endl;

             int sum=0;
             for(int z=idx_info; z>=0; z--){
                 cout<<bit_rep[z];
                 sum=sum+bit_rep[z]*pow(2,z);
             }
             cout<<endl<<"sum: "<<sum<<endl;

             break;
         }
         idx_info--;
     }

     idx_info=31;

     int m=ceil(log2(idx_info+1));
     int n=m+idx_info+1;

     int parity[10]; // stores the locations of the parity bits

     // generate the position of the parity bits
     parity[0] = 1;
     for(int i=1; i<10; i++){
         parity[i]=(parity[i-1]<<1) & 0xfffffffe;
     }

    // compute the value of the bits at those positions
    int red[1024];
    int l=0;
    int result=0;
    int test=0;
    for (int j=0; j<m; j++)
    {
      red[j] = 0;
      l = 0;
      for (int i=0; i<n; i++)
      {
        // Check that "i" is not a parity position = not a power of 2
        result = 0;
        test = 1;
        for (int idx=0; idx<m; idx++)
        {
          if (i==test) result = 1;
          test *= 2;
        }
        if (!result)
        {
          l++;
          if ( (i>>(j)) & 0x01 )
            red[j] ^= bit_rep[l];
        }
      }
    }
    printf("parity positions: ");
    for (int i=0; i<m; i++) printf("%2d ", parity[i]);
    printf("\n");

    printf("parity bits = ");
    for (int j=0; j<m; j++) printf("%1d", red[j]);
    printf("\n");

     // Transmit codeword
     int i=0;
     l=0;

     vector<int> code(n, 0);
     for(int j=0; j<n; j++){
        if(j==(parity[l]-1) && l<m){
            code[j]=red[l];
            l++;
        }
        else{
            code[j]=bit_rep[idx_info-i];
            i++;
        }
    }


     return code;
 }

 /* decode
  * DESCRIPTION: Decodes a binary hamming code output
  * INPUT: A binary encoded hamming code thing
  * OUTPUT: The original number which was encoded using the
  */
  double decode(vector<int> encoded){
    // compute syndrome
    int syn = 0;
    int n=encoded.size()-1;
    for(int i=1; i<=n; i++)
        if (encoded[i])
            syn ^= i;

    // correct error if needed
    if(syn)
        encoded[syn]^=1;

    // compute estimate
    int sum=0;
    cout<<n<<endl;
    for(int z=1; z<=n; z++){
        // cout<<encoded[z];
        if((int)log2(z)!=log2(z)){
            cout<<encoded[z]<<" "<<z<<endl;
        }
    }

    cout<<"estimated sum: "<<sum<<endl;

    return 0;
  }

/* R
 * DESCRIPTION: This is a basic LDP randomizer this is based on alogrithm 1 of
 * the paper
 * INPUTS: A vector of vertices on the hypercube or an empty vector of zeros
 * NOTE: the input values are always going to be {0, 1/sqrt(x.size()), -1/sqrt(x.size())}
 * OUTPUTS: vector double same size as input x but with only one nonzero vlaue
 */
pair<int, double> R(std::vector<int> x, double epsilon){
    if(x.size()<=0 || epsilon<=0){
        return pair<int, double>(-1,-1);
    }

    // BEGIN ALGORITHM
    int m=m_code; //x.size();
    // time_t timer; // TODO: figure how to make this random
    default_random_engine gen(time(NULL));

    //////////////////////////////  STEP 1
    // sample j from m uniformally at random where j is an index into
    // the vector m thus having index value of 0 to m.size()-1
    uniform_int_distribution<int> idx_gen(0, m-1);
    int j=idx_gen(gen);

    //////////////////////////////  STEP 2
    // if x[j]!=0
    // randomize jth input to value to {-1/sqrt(m), 1/sqrt(m)}
    // with the following distribution
    //               exp(epsilon)+1
    // let c_epsilon=--------------= O(1/epsilon)
    //               exp(epsilon)-1
    //      --
    //      | c_epsilon * m * x_j,    w.p. exp(epsilon)/(exp(epsilon)+1)
    // zj = |
    //      | -c_epsilon * m * x_j,    w.p. 1/(exp(epsilon)+1)
    //      --
    // NOTE: x_j is going to be +- 1/sqrt(m) therefore it's value is still going
    // to be the same as the else case
    double zj;
    double c_epsilon=(exp(epsilon)+1)/(exp(epsilon)-1);

    if(x[j]!=0){
        bernoulli_distribution zj_chooser(exp(epsilon)/(exp(epsilon)+1));

        if(zj_chooser(gen)){// positive choice
            zj=c_epsilon*m*x[j];
        }
        else{ // negative choice
            zj=-c_epsilon*m*x[j];
        }
    }
    // else generate uniform bit z_j={c_epsilon*sqrt(m), -c_epsilon*sqrt(m)}
    else{
        uniform_int_distribution<int> _0_1(0,1);
        int res=_0_1(gen);

        if(res==0){
            zj=c_epsilon*sqrt(m);
        }
        else{
            zj=c_epsilon*sqrt(m)*-1.0;
        }
    }

    //////////////////////////////  STEP 3
    // return only the pair of index and the value
    pair<int, double> result; // first=idx, second=value
    result.first=j;
    result.second=zj;
    vector<double> z_(m, 0.0);

    return result;
}

/* GenProj
 * DESCRIPTION: Generates a random projection matrix following the
 * Johnson-Lindenstrauss lemma
 * INPUT: Two integers representing the size m=rows, d=cols
 * OUTPUT: A Projection matrix with sufficient randomness? Something...
 */
 dataSet GenProj(int m, int d){
     dataSet result;

     // Uniformally generate a binary value either 0=1/sqrt(m) and 1=-1/sqrt(m)
     // Then proceed to assign the value to each thing in the matrix

     return result;
 }
/* PROT_FO
 * DESCRIPTION: This creates a frequency oracle from a list of user inputs
 * and a confidence parameter
 * INPUT: Input vector of user data, and a confidence parameter
 * OUTPUT: A dataSet structure that contains the phi FO matrix
 */
dataSet PROT_FO(std::vector<double> v, double epsilon, double beta){
    dataSet FO;

    /////////////////// STEP 1 Initialize the values required for computation
    int d; // TODO: Need to figure out where to set this this is the universe
           // of possible values (aka number of values possible)
    int n=v.size();
    double gamma=sqrt(log(2*d/beta)/(epsilon*n));
    int m=(log(d+1)*log(2/beta))/(gamma*gamma);

    ////////////////// STEP 2 Create a phi matrix
    dataSet phi=GenProj(m,d);

    return FO;
}

/* PROT_PP_S_Hist_PP
 * DESCRIPTION: Produces a succinct histogram under promise
 * INPUTS: data each row=1 user each column = attribute, epsilon (privacy),
 * and beta (confidence parameter)
 * OUTPUTS:
 */
pair<double, double> PROT_PP_S_Hist_pp(dataSet &data, double epsilon){
    /////////////////////////// STEP 1
    // storage for the value and index into the {0,0...,zj,...0,0} vector
    vector<pair<int, double>> ind_val;
    // each user will encode their item
    for(int i=0; i<data.get_m(); i++){
        // if item exists code
        vector<int> encoded(m_code, 0);
        if(data.get_dataMtx(i, 0)){
            encoded=code((int)data.get_dataMtx(i, 0));
        }
        // else set zero
        else{}
        // compute user randomized report using R and sends to server
        ind_val.push_back(R(encoded, epsilon));
    }

    /////////////////////////// STEP 2
    // server computes z_bar average over all the z's
    vector<double> sums(m_code, 0.0);
    vector<double> counts(m_code, 0);
    for(int i=0; i<(int)ind_val.size(); i++){
        sums[ind_val[i].first]+=sums[ind_val[i].first]+ind_val[i].second;
        counts[ind_val[i].first]++;
    }
    vector<double> avgs(m_code, 0.0);  // this is z_bar
    for(int i=0; i<(int)sums.size(); i++){
        avgs[i]=(double)sums[i]/counts[i];
    }

    // compute y by rounding to the 1/sqrt(m) thing
    vector<double> y_bar(m_code, 0.0);
    for(int i=0; i<(int)avgs.size(); i++){
        if(y_bar[i]>=0)
            y_bar[i]=1.0/sqrt(m_code);
        else
            y_bar[i]=-1.0/sqrt(m_code);
    }

    /////////////////////////// STEP 3
    // decode the y vector to get the average and the frequency estimate
    // convert y_bar to binary
    vector<int> y_bar_int(m_code, 0);
    for(int i=0; i<(int)y_bar.size(); i++){
        if(y_bar[i]==1.0/sqrt(m_code))
            y_bar_int[i]=1;
        else
            y_bar_int[i]=0;
    }
    double avg_est=decode(y_bar_int);

    int inner_product=0.0; // this is the frequency estimate
    for(int i=0; i<(int)avgs.size(); i++){
        inner_product+=(double)avgs[i]*y_bar[i];
    }

    /////////////////////////// STEP 5
    // return the estimated hevay hitter value and its frequency
    pair<double, double> result;
    result.first=avg_est;
    result.second=inner_product;

    return result;
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
