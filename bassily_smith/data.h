#ifndef __DATA_H
#define __DATA_H

#include <vector>

class boolSet;

class dataSet{
    public:
        dataSet(); // constructor

        // data matrix getter/setter
        double get_dataMtx(int m_in, int n_in);
        int set_dataMtx(int m_in, int n_in, double val);

        // matrix attributes getter
        int get_m();
        int get_n();

        // resizing the data matrix
        int rekeep_dataMtx(int m_in, int n_in);

        // resize the data matrix and keep data
        int resize_dataMtx(int m_in, int n_in);

        // edit column/attribute
        int columnSet(int n_in, std::vector<double> perturbedColumn);


    private:
        int m, n;
        std::vector< std::vector<double> > dataMtx;
};

class boolSet{
    public:
        boolSet(); // constructor

        // data matrix getter/setter
        bool get_boolMtx(int m_in, int n_in);
        int set_boolMtx(int m_in, int n_in, bool val);

        // matrix attributes getter
        int get_m();
        int get_n();

        // resizing the data matrix
        int rekeep_boolMtx(int m_in, int n_in);

        // resize the data matrix and keep data
        int resize_boolMtx(int m_in, int n_in);

        // edit column/attribute
        int columnSet(int n_in, std::vector<double> perturbedColumn);


    private:
        int m, n;
        std::vector< std::vector<bool> > dataMtx;
};

#endif
