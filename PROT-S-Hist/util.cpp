#include "util.h"
#include <random>
#include <iomanip> // remove later
#include <cstdlib>
#include <ctime>
#include <bitset>
#include "c++/ezpwd/bch"

#define TAB 0x09
#define LARGENUMBER 9999999999
int m_max=7; // set number of bits of R
using namespace std;
static default_random_engine gen;

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
 * DESCRIPTION: Encodes a thing into binary hamming code (7,4)
 * INPUT: A number to encode into binary hamming code
 * OUTPUT: A character array of the encoded hamming stuff
 */
 vector<int> code(int input){
     // deal with 0 input case
     if(input==0){
         vector<int> asdf(m_max+3, 0);
         return asdf;
    }
     // create binary input, figure out the length of it
     // create hamming code as required

     bitset<32> bit_rep(input); // displays input as binary representation
                             // note the bit_rep[0] repreents the smallest ibt
                             // this means that 11 -> 1011 in binary
                             // bit_rep[3]=1; bit_rep[2]=0; bit_rep[1]=1 bit_rep[0]=1;
     int right, left; // 01100 read right to left so init, right=1 left=0
                      // then right=1 left=1 ---->
     int idx_info=31;
     while(1){
         right=bit_rep[idx_info-1];
         left=bit_rep[idx_info];

         if(right==1 && left==0){
             //cout<<"input: "<<input<<endl;
             idx_info=idx_info-1; // this is now the number of bits
             //cout<<"number of digits: "<<idx_info+1<<endl;

             int sum=0;
             for(int z=idx_info; z>=0; z--){
                 //cout<<bit_rep[z];
                 sum=sum+bit_rep[z]*pow(2,z);
             }
             //cout<<endl<<"sum: "<<sum<<endl;

             break;
         }
         idx_info--;
     }

     // we now create 32/4=8 sets of 7 bits
     // number of sets of 7 we need
     int num_blocks=idx_info/4+1; int block_counter=num_blocks;
     vector<int> hamming_code(num_blocks*7, 0); // store the bits backwards
                                    // this means that bit_rep[31]==hamming_code[31]

    // the following computation is computed by doing
    //                    p1 p2 p3 d1 d2 d3 d4 vectors
    //                   | 0  1  1  1  0  0  0 |
    // | d1 d2 d3 d4 | * | 1  0  1  0  1  0  0 | = encoded
    //                   | 1  1  0  0  0  1  0 |
    //                   | 1  1  1  0  0  0  1 |
     for(int i=0; i<num_blocks; i++){
        // 0 1 1 1
        hamming_code[7*i]=(bit_rep[4*block_counter-1]*0+bit_rep[4*block_counter-2]*1
                         +bit_rep[4*block_counter-3]*1+bit_rep[4*block_counter-4]*1)%2;

        // 1 0 1 1
        hamming_code[7*i+1]=(bit_rep[4*block_counter-1]*1+bit_rep[4*block_counter-2]*0
                         +bit_rep[4*block_counter-3]*1+bit_rep[4*block_counter-4]*1)%2;

        // 1 1 0 1
        hamming_code[7*i+2]=(bit_rep[4*block_counter-1]*1+bit_rep[4*block_counter-2]*1
                         +bit_rep[4*block_counter-3]*0+bit_rep[4*block_counter-4]*1)%2;

        // 1 0 0 0
        hamming_code[7*i+3]=(bit_rep[4*block_counter-1]*1+bit_rep[4*block_counter-2]*0
                         +bit_rep[4*block_counter-3]*0+bit_rep[4*block_counter-4]*0)%2;

        // 0 1 0 0
        hamming_code[7*i+4]=(bit_rep[4*block_counter-1]*0+bit_rep[4*block_counter-2]*1
                         +bit_rep[4*block_counter-3]*0+bit_rep[4*block_counter-4]*0)%2;

        // 0 0 1 0
        hamming_code[7*i+5]=(bit_rep[4*block_counter-1]*0+bit_rep[4*block_counter-2]*0
                         +bit_rep[4*block_counter-3]*1+bit_rep[4*block_counter-4]*0)%2;

        // 0 0 0 1
        hamming_code[7*i+6]=(bit_rep[4*block_counter-1]*0+bit_rep[4*block_counter-2]*0
                         +bit_rep[4*block_counter-3]*0+bit_rep[4*block_counter-4]*1)%2;

     }

     return hamming_code;
 }

 /* decode
  * DESCRIPTION: Decodes a binary hamming code output
  * INPUT: A binary encoded hamming code thing
  * OUTPUT: The original number which was encoded using the
  */
  double decode(vector<int> encoded){
      // check parity bits to see if they match
      int counter=encoded.size()/7; // used to mark how many parts of 7

      // create a parity check matrix and solve for the syndrome
      // Row 1: Contains 1 in the first parity bit position and the bits used to calculate parity
      // Row 2: Contains 1 in the second parity bit position and the bits used to calculate parity
      // Row 3: Contains 1 in the third partiy bit position and the bits used to calculate parity
      //
      //     | 1 0 0 0 1 1 1 |
      // H = | 0 1 0 1 0 1 1 |
      //     | 0 0 1 1 1 0 1 |
      //
      // H*encode= 3x1 matrix of syndrome
      // if syndrome all 0 then error free
      // else flipping the encoded bit that is in the position of the column in [H] that matches the syndrome will result in a valid code word
      dataSet *H=new(dataSet);
      H->rekeep_dataMtx(3, 7);
      // Row 1
      H->set_dataMtx(0+1, 0+1, 1); H->set_dataMtx(0+1, 1+1, 0); H->set_dataMtx(0+1, 2+1, 0);
      H->set_dataMtx(0+1, 3+1, 0); H->set_dataMtx(0+1, 4+1, 1); H->set_dataMtx(0+1, 5+1, 1);
      H->set_dataMtx(0+1, 6+1, 1);
      // Row 2
      H->set_dataMtx(1+1, 0+1, 0); H->set_dataMtx(1+1, 1+1, 1); H->set_dataMtx(1+1, 2+1, 0);
      H->set_dataMtx(1+1, 3+1, 1); H->set_dataMtx(1+1, 4+1, 0); H->set_dataMtx(1+1, 5+1, 1);
      H->set_dataMtx(1+1, 6+1, 1);
      // Row 3
      H->set_dataMtx(2+1, 0+1, 0); H->set_dataMtx(2+1, 1+1, 0); H->set_dataMtx(2+1, 2+1, 1);
      H->set_dataMtx(2+1, 3+1, 1); H->set_dataMtx(2+1, 4+1, 1); H->set_dataMtx(2+1, 5+1, 0);
      H->set_dataMtx(2+1, 6+1, 1);

      // Multiply H with encoded word
      dataSet syndrome;
      syndrome.rekeep_dataMtx(3, 1); // one row for each bit
      double sum=0.0; // used to keep track of the sum we need to return the decoded value as
      int pow_count=0; // used to count the powers of two we're adding together
      for(int wc=0; wc<(int)encoded.size()/7; wc++){
          // H*encode is to create the syndrome
          // 3 x 7 * 7 x 1 -> 3x1
          int temp;
          temp=H->get_dataMtx(0, 0)*encoded[wc*7+0]+H->get_dataMtx(0, 1)*encoded[wc*7+1]
            +H->get_dataMtx(0, 2)*encoded[wc*7+2]+H->get_dataMtx(0, 3)*encoded[wc*7+3]
            +H->get_dataMtx(0, 4)*encoded[wc*7+4]+H->get_dataMtx(0, 5)*encoded[wc*7+5]
            +H->get_dataMtx(0, 6)*encoded[wc*7+6];
          temp=temp%2;
          syndrome.set_dataMtx(0+1, 0+1, temp);

          temp=H->get_dataMtx(1, 0)*encoded[wc*7+0]+H->get_dataMtx(1, 1)*encoded[wc*7+1]
            +H->get_dataMtx(1, 2)*encoded[wc*7+2]+H->get_dataMtx(1, 3)*encoded[wc*7+3]
            +H->get_dataMtx(1, 4)*encoded[wc*7+4]+H->get_dataMtx(1, 5)*encoded[wc*7+5]
            +H->get_dataMtx(1, 6)*encoded[wc*7+6];
          temp=temp%2;
          syndrome.set_dataMtx(1+1, 0+1, temp);

          temp=H->get_dataMtx(2, 0)*encoded[wc*7+0]+H->get_dataMtx(2, 1)*encoded[wc*7+1]
            +H->get_dataMtx(2, 2)*encoded[wc*7+2]+H->get_dataMtx(2, 3)*encoded[wc*7+3]
            +H->get_dataMtx(2, 4)*encoded[wc*7+4]+H->get_dataMtx(2, 5)*encoded[wc*7+5]
            +H->get_dataMtx(2, 6)*encoded[wc*7+6];
          temp=temp%2;
          syndrome.set_dataMtx(2+1, 0+1, temp);

          // do we have errors? if so correct
          if(syndrome.get_dataMtx(0, 0)!=0 || syndrome.get_dataMtx(1, 0)!=0 || syndrome.get_dataMtx(2, 0)!=0){
              // if the columns dont sum to zero then we need to flip that bit
              int bit_col;
              for(int a=0; a<7; a++){
                  int col_sum=0;
                  col_sum+=H->get_dataMtx(0, a)*encoded[wc*7+a]+H->get_dataMtx(1, a)*encoded[wc*7+a]
                    +H->get_dataMtx(2, a)*encoded[wc*7+a];
                  if(col_sum>0){
                      bit_col=a; // this is now the column we need to fix
                  }
              }
              if(encoded[wc*7+bit_col]==0){encoded[wc*7+bit_col]=1;}
              else{encoded[wc*7+bit_col]=0;}
          }
          // else proceed
          for(int p=3; p>=0; p--){
              sum+=pow(2, pow_count)*encoded[wc*7+3+p]; // word_offset+parity_offset+offset
              pow_count++;
          }
      }

    //   cout<<"post syndrome"<<endl;
    //   for(int i=0; i<(int)encoded.size(); i++){
    //       cout<<encoded[i];
    //   }
    //   cout<<endl;

      delete H;
      return sum;
  }

/* R
 * DESCRIPTION: This is a basic LDP randomizer this is based on alogrithm 1 of
 * the paper
 * INPUTS: A vector of vertices on the hypercube or an empty vector of zeros
 * NOTE: the input values are always going to be {0, 1/sqrt(x.size()), -1/sqrt(x.size())}
 * OUTPUTS: vector double same size as input x but with only one nonzero vlaue
 */
pair<int, double> R(std::vector<int> &x, double epsilon){
    if(x.size()<=0 || epsilon<=0){
        return pair<int, double>(-1,-1);
    }

    // BEGIN ALGORITHM
    int m=m_max-1; //x.size();
    // time_t timer;
    // static default_random_engine gen;//(time(0));
    // static default_random_engine gen2;

    // VERIFY THAT THE VECTOR IS COMPLETELY EMPTY OR NOT
    int flag=0; // if flag==1 then we have non zero values
    for(int i=0; i<(int)x.size(); i++){
        if(x[i]!=0){
            flag=1;
        }
    }

    //////////////////////////////  STEP 1
    // sample j from m uniformally at random where j is an index into
    // the vector m thus having index value of 0 to m.size()-1
    // uniform_int_distribution<int> idx_gen(0, m);
    int j=rand()%m; // idx_gen(gen);

    //////////////////////////////  STEP 2
    // if x=/=0
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

    if(flag){
        bernoulli_distribution zj_chooser(exp(epsilon)/(exp(epsilon)+1));

        if(zj_chooser(gen)){// positive choice
            zj=c_epsilon*sqrt(m+1);//c_epsilon*(m+1)*x[j];
        }
        else{ // negative choice
            zj=-c_epsilon*sqrt(m+1);//-c_epsilon*(m+1)*x[j];
        }
    }
    // else generate uniform bit z_j={c_epsilon*sqrt(m), -c_epsilon*sqrt(m)}
    else{
        uniform_int_distribution<int> _0_1(0,1);
        int res=_0_1(gen);

        if(res==0){
            zj=c_epsilon*sqrt(m+1);
        }
        else{
            zj=c_epsilon*sqrt(m+1)*-1.0;
        }
    }

    //////////////////////////////  STEP 3
    // return only the pair of index and the value
    pair<int, double> result; // first=idx, second=value
    result.first=j;
    result.second=zj;
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
    pair<double, double> result; // return values, first is the value
                                // second is the frequency
    int _0=0; int _1=0; int total=0;
    ///////////////////////// STEP 1
    // 1: for Users i = 1 to n do
    // 2:   If vi =/= NOT, then user i encodes its item: xi = c(vi). Else, user i sets xi = 0.
    // 3:   User i computes its private report: zi = Ri (xi, epsilon).
    // 4:   User i sends zi the server
    vector<pair<int, double> > zi_temp;
    pair<int,double> zi;
    for(int i=0; i<data.get_m(); i++){
        if(data.get_dataMtx(i, 0)!=0){ // code
            vector<int> c=code(data.get_dataMtx(i,0));
            //if((int)c.size()>m_max){m_max=c.size();}
            //cout<<c.size()<<" "<<m_max<<endl;
            zi=R(c, epsilon);
            if(zi.second>0){_0++;total++;}
            else{_1++;total++;}
            zi_temp.push_back(zi);
        }
        else{ // set 0
            vector<int> c(7,0);
            zi=R(c, epsilon);
            if(zi.second>0){_0++;total++;}
            else{_1++;total++;}
            zi_temp.push_back(zi);
        }
        //cout<<"i: "<<i<<" "<<"zi (idx, val): "<<zi.first<<" "<<zi.second<<endl;
    }
    cout<<"0%: "<<(double)_0/total<<" total: "<<total<<" 0: "<<_0<<endl;

    ///////////////////////// STEP 2
    // 5: Server computes z_bar=(1/n)summation(1 to n)zi
    // now actually store all the data in a zi vector
    // zi.first is idx, zi.second is value
    // zi is stored in the zi_temp vector
    vector<double> z(m_max, 0.0);
    for(int i=0; i<(int)zi_temp.size(); i++){
        z[zi_temp[i].first]=z[zi_temp[i].first]+zi_temp[i].second;
    }
    for(int i=0; i<(int)m_max; i++){
        z[i]=z[i]/zi_temp.size(); // this should be == zi_temp.size()
    }

    ///////////////////////// STEP 3
    // compute y by rounding z to {-1/sqrt(m), 1/sqrt(m)}^m
    //       --
    //       | 1/sqrt(m)    if zj>=0
    // y_j = |
    //       | -1/sqrt(m)   otherwise
    //       --
    vector<int> y(z.size(), 0.0);
    for(int i=0; i<(int)z.size(); i++){
        if(z[i]>0)
            y[i]=1;//1/sqrt(m_max);
        else
            y[i]=0;//-1/sqrt(m_max);
    }

    ///////////////////////// STEP 4
    // decode y
    result.first=decode(y);

    vector<double> y_d(z.size(), 0.0);
    vector<int> v_hat=code(result.first);
    vector<double> v_hat_d(z.size(), 0.0);
    for(int i=0; i<(int)z.size(); i++){
        if(z[i]>0)
            y_d[i]=1/sqrt(m_max);
        else if(z[i]<=0)
            y_d[i]=-1/sqrt(m_max);
        if(v_hat[i]==1)
            v_hat_d[i]=1/sqrt(m_max);
        else if(v_hat[i]==0)
            v_hat_d[i]=-1.0/sqrt(m_max);

    }

    double freq=0.0;
    for(int i=0; i<(int)v_hat.size(); i++){
        freq+=v_hat_d[i]*y_d[i];
    }
    result.second=freq;

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
