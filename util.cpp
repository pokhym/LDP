#include "util.h"

#define TAB 0x09

using namespace std;

/* parseData
 * DESCRIPTION: Parses data read from a txt file exported from a CSV in excel
 * assumes that the columns are tab separated
 * INPUTS:  First column=names of itmes
 *          All columns thereafter should contain some form of number
 * OUTPUTS: 0 on success and -1 on fail
 */
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
                            //cout<<"row_count: "<<row_count<<endl;
                            data.rekeep_dataMtx(row_count, max_col);
                            data.set_dataMtx(row_count, col_count, 0.0);
                            data.set_dataMtx(0, row_count, row_count);
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
                                    
                                    cout<<"before dec after dec"<<before_decimal<<" "<<after_decimal<<endl;

                                    // calculate value before decimal
                                    for(int count4=0; count4<before_decimal; count4++){
                                        val=val+(double)pow(10,count4)*(temp_num_str[count4]-0x30);
                                    }

                                    // calculate value after decimal
                                    for(int count4=1; count4<=after_decimal; count4++){
                                       val=val+(double)(1.0/pow(10,count4))*
                                           1.0*(temp_num_str[count4+before_decimal]-0x30);

                                    }

                                    cout<<val<<endl;

                                    // save into data matrix
                                    //cout<<"row_count: "<<row_count<<" "<<col_count<<endl;
                                    data.rekeep_dataMtx(row_count, max_col);
                                    data.set_dataMtx(row_count, col_count, val);
                                    data.set_dataMtx(0, row_count, row_count);

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
                cout<<temp_str[j];
                j++;
            }
            if(col_count>max_col)
                max_col=col_count;
            col_count=0; // reset column count
        }
    }
    // read the number of rows and columns
    // update m, n in the data set
    // recall the column=attribute and row=value
    return 0;
}
