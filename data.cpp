#include <iostream>
#include <vector>
#include "data.h"

/* dataSet()
 * DESCRIPTION: Default constructor
 * INPUTS: None
 * OUTPUTS: None
 */
dataSet::dataSet(){
    m=1;
    n=1;

    // resize col, row
    dataMtx.resize(n, std::vector<double>(m));
    dataMtx[0][0]=1.0;
}

/* get_dataMtx
 * DESCRIPTION: returns the matrix value at an index
 * INPUTS: indices
 * OUTPUTS: matrix value at index or -1 on fail
 */
double dataSet::get_dataMtx(int m_in, int n_in){
    // error checkgin
    if(m_in<0 || m_in>m || n_in<0 || n_in>n)
        return -1;

    return dataMtx[m_in][n_in];
}

/* set_dataMtx
 * DESCRIPTION: sets a index in the matrix
 * INPUTS: m_in, n_in two indices
 * OUTPUTS: 0 success -1 faile
 */
int dataSet::set_dataMtx(int m_in, int n_in, double val){
    // check bounds
    if(m_in<0 || m_in>this->m || n_in<0 || n_in>this->n)
        return -1;

    dataMtx[m_in-1][n_in-1]=val;
    return 0;
}

/* get_m
 * DESCRIPTION: gets the m value
 * INPUTS: none
 * OUTPUTS m value
 */
int dataSet::get_m(){
    return m;
}

/* get_n
 * DESCRIPTION: gets the n value
 * INPUTS: none
 * OUTPUTS: n value
 */
int dataSet::get_n(){
    return n;
}

/* rekeep_dataMtx
 * DESCRIPTION: this function resizes the matrix while keeping the previous data intact
 * INPUTS: the new size of the matrix
 * OUTPUTS: 0 success -1 fail
 */
int dataSet::rekeep_dataMtx(int m_in, int n_in){
    // check bounds
    if(m_in==m && n_in==n)
        return 0;
    else if(m_in<0 || n_in<0 || m_in<m || n_in<n)
        return -1;

    // resize automatically retains data
    // v seems broken it doesnt seem to take into account dynamic allocation need to do manual crap
    // dataMtx.resize(m_in, std::vector<double>(n_in));
 
    // generate new vector
    if(m_in>m || n_in>n){
        std::vector< std::vector<double> > asdf(m_in*2, std::vector<double>(n_in*2, 0.0) );
    
        // copy into new one and reassign
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                asdf[i][j]=dataMtx[i][j];
            }
        }
        dataMtx=asdf;

        m=m_in;
        n=n_in;

        return 0;
    }
    else{
        return 0;
    }

}

/* resize_dataMtx
 * DESCRIPTION: this function resizes and clears the matrix
 * INPUTS: new size of matrix
 * OUTPUTS: 0 on success -1 on fail
 */
int dataSet::resize_dataMtx(int m_in, int n_in){
    // check bounds
    if(m_in==m && n_in==n){
        // extra clear here
        for(int i=0; i<m_in; i++){
            for(int j=0; j<n_in; j++){
                dataMtx[i][j]=0.0;
            }
        }
        return 0;
    }
    else if(m_in<0 || n_in<0)
        return -1;
    
    // resize
   dataMtx.resize(n_in, std::vector<double>(m_in));
    // clear
     for(int i=0; i<m_in; i++){
        for(int j=0; j<n_in; j++){
           dataMtx[i][j]=0.0;
        }
    }
     m=m_in;
     n=n_in;
    return 0;
}

/* columnEdit
 * DESCRIPTION: Edits one column of the dataMtx
 * INPUTS: column to edit, new data
 * OUTPUTS: 0 on success -1 on fail
 */
int dataSet::columnSet(int n_in, std::vector<double> perturbedColumn){
    n_in=n_in+1;
    if(perturbedColumn.size()!=(unsigned int)m|| n_in<=0 || n_in>n)
        return -1;
    
    // loop through all the values
    for(int i=0; i<m; i++){
        this->set_dataMtx(i+1, n_in, perturbedColumn[i]);
    }
    
    return 0;
}
