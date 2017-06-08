#include "util.h"
#include <random>

#define TAB 0x09
#define LARGENUMBER 9999999999
#define CR 0xD

using namespace std;

/* parseData
 * DESCRIPTION: Parses a comma separated list of data. Each column represents
 * an attribute an each row represents a user
 * INPUTS: filename, dataset
 * OUTPUTS: 0 on success -1 on failure
 */
int parseData(char const* filename, dataSet &data){
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
    
    for(int i=0; i<data.get_m(); i++){
        for(int j=0; j<data.get_n(); j++){
            double asdf=data.get_dataMtx(i,j);
            cout<<asdf<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    
    myfile.close();
    
    return 0;
}

/* normalizeNeg1toPos1
 * DESCRIPTION: normalizes each column of the input data set -1 to 1
 * INPUTS: dataSet, outlier condition
 * OUTPUTS: 0 on success -1 on fail
 */
int normalizeNeg1toPos1(dataSet &data, vector<double> outlier){
    if(data.get_m()==1 || data.get_n()==1 || outlier.size()==0)
        return -1;
    
    // store the max values of each for normalization
    vector<double> max(data.get_n(), -LARGENUMBER);
    vector<double> min(data.get_n(), LARGENUMBER);

    // loop through rows and columns to obtain max
    double max_curr_col=-LARGENUMBER; // maximum in current column
    double min_curr_col=LARGENUMBER; // minimum in current column
    for(int n=1; n<data.get_n(); n++){
        for(int m=0; m<data.get_m(); m++){
            if(data.get_dataMtx(m,n)!=outlier[n]){
                if(data.get_dataMtx(m,n)>max_curr_col)
                    max_curr_col=data.get_dataMtx(m,n); // update

                if(data.get_dataMtx(m,n)<min_curr_col)
                    min_curr_col=data.get_dataMtx(m,n); // update
            }
        }
        // save into vector and iterate to next column
        max[n]=max_curr_col;
        min[n]=min_curr_col;
    }
    
    // normalize per column
    for(int n=1; n<data.get_n(); n++){
        for(int m=0; m<data.get_m(); m++){
            // if maximum is 0 (all outliers) don't do anything
            if(max[n]==0){}

            // else we need to normalize all nonzero values
            else{
                // if the data in that m,n is not an outlier normalize
                if(data.get_dataMtx(m,n)!=outlier[n]){
                    double new_val=data.get_dataMtx(m,n);
                    new_val=(new_val-min[n])/(max[n]-min[n]);
                    
                    // center around -1 to 1
                    new_val=new_val-0.5;
                    new_val=new_val*2.0;

                    // save new value into data
                    data.set_dataMtx(m,n,new_val);
                }

                // else we need to set to 0
                else{
                    data.set_dataMtx(m,n,0.0);
                }
            }
        }
    }
    return 0;
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
    int n=2;
    
    // loop through all the different entries (not attributes)
    for(int i=0; i<data.get_m(); i++){
        double probability, numerator, denominator;
        //          ti[Aj]*(exp^(epsilon)-1)+exp^(epsilon)+1
        // Pr[u=1]= ----------------------------------------
        //                     2*exp^(epsilon)+2
        
        numerator=data.get_dataMtx(i, n)*(exp(epsilon)-1)+exp(epsilon)+1;
        denominator=(double)exp(epsilon)+2;
        probability=numerator/denominator;
        
        // generate random crap code pulled crom c++reference
        random_device rd;
        mt19937 gen(rd());
        bernoulli_distribution d(probability); // set p

        // grab data values from bernoulli
        int bernoulli_result=d(gen);
        
        // assign values
        if(bernoulli_result){
            perturbed[i]=(double)(data.get_n()-1)*(exp(epsilon)+1.0)/(exp(epsilon)-1.0);
        }
        else{
            perturbed[i]=(double)-1.0*(data.get_n()-1)*(exp(epsilon)+1.0)/(exp(epsilon)-1.0);
        }
    }

    // update matrix with the new perturbed column
    data.columnSet(n, perturbed);
    
    return perturbed;
}

/* parseData
 * DESCRIPTION: Parses data read from a txt file exported from a CSV in excel
 * assumes that the columns are tab separated
 * INPUTS:  First column=names of itmes
 *          All columns thereafter should contain some form of number
 * OUTPUTS: 0 on success and -1 on fail
 */

/*
int parseData(char const* filename, dataSet &data){
    string line;
    ifstream myfile(filename);
    int col_count=0;
    int max_col=0;
    int row_count=0;

    // open file
    if(myfile.is_open()){

        // read all lines
        while(getline(myfile, line)){
            // increment row, column counts
            row_count++;

            // obtain length of line
            int str_len=line.size();
            string temp_str;

            // copy the string into a the temp until hit tab
            int j=0;
            for(int i=0; i<str_len; i++){
                
                // if we hit a tab process and clear string
                if(line[i]==TAB || i==str_len-1){
                    col_count++;
                    if(col_count>max_col)
                        max_col=col_count;
                    // case where we have parsed the name
                    if(col_count==1){
                        j=0; i++;
                        temp_str.clear();
                    }

                    // else we have a data column
                    else{
                        if(temp_str[0]=='F'){
                            // save into data matrix
                            data.rekeep_dataMtx(row_count, max_col);
                            data.set_dataMtx(row_count-1, col_count-1, 0.0);
                            data.set_dataMtx(row_count-1, 0, row_count);
                        }
                        else{
                            for(int count=0; count<j; count++){                            
                                // take everything after the dollar sign (case
                                // price)
                                if(temp_str[count]=='$'){
                                    string temp_num_str;
                                    int count2;
                                    int count3=0;

                                    // grab only the number part
                                    for(count2=(count+1); count2<j; count2++){
                                        temp_num_str[count3]=temp_str[count2];
                                        count3++;
                                    }

                                    // convert number part into decimal
                                    double val;
                                    int before_decimal, after_decimal, decimal_flag;
                                    before_decimal=0; after_decimal=0; decimal_flag=0;
                                    
                                    for(int count4=0; count4<count3; count4++){
                                        if(temp_num_str[count4]=='.'){
                                            count4++;
                                            decimal_flag=1;
                                        }
                                        if(temp_num_str[count4]==' ') // weird error case
                                            break;
                                        if(decimal_flag==0)
                                            before_decimal++;
                                        if(decimal_flag==1)
                                            after_decimal++;
                                    }
                                    
                                    // calculate value before decimal
                                    for(int count4=0; count4<before_decimal; count4++){
                                        val=val+(double)pow(10,count4)*(temp_num_str[count4]-0x30);
                                    }

                                    // calculate value after decimal
                                    for(int count4=1; count4<=after_decimal; count4++){
                                       val=val+(double)(1.0/pow(10,count4))*
                                           1.0*(temp_num_str[count4+before_decimal]-0x30);

                                    }

                                    // save into data matrix
                                    data.rekeep_dataMtx(row_count, max_col);
                                    data.set_dataMtx(row_count-1, col_count-1, val);
                                    data.set_dataMtx(row_count-1, 0, row_count);

                                    // reset values
                                    val=0;
                                    j=0;
                                    temp_str.clear();
                                }
                            }
                        }
                    }
                }
                temp_str[j]=line[i]; // copy string
                j++;
            }
            if(col_count>max_col)
                max_col=col_count;
            col_count=0; // reset column count
        }

    myfile.close();
    }
    // read the number of rows and columns
    // update m, n in the data set
    // recall the column=attribute and row=value
    return 0;
}
*/
