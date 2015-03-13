#ifndef SPHERICITY_h
#define SPHERICITY_h

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <stdexcept>
#include <exception>
#include <TROOT.h>
#include <vector>
#include <sstream>
#include <TString.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TArrayD.h>
#include <TVectorD.h>
using namespace std;


//////////////////////////////////////////////////////
//
//this fucntion is to calculate thrust and transverse thrust ////
//please take care that this function returns minus of thrust
//because of Minimizer
/////////////////////////////////////////////////////

class Sphericity{

public:
    Sphericity(){Initialize();};
    ~Sphericity(){};
    ////for start sphericity tool///
    void Initialize(){ m_v_tlv.clear(); m_v_tlv_size=0;};
    void SetTLV(const vector<TLorentzVector> &v_tlv, int v_tlv_size);
    ////for calculate and get the result//
    
    void GetSphericity(Double_t &S, Double_t &ST, Double_t &A);
    void GetBoostedSphericity(Double_t &S, double &ST, Double_t &A);
    void GetSafeSphericity(Double_t &S, Double_t &ST, Double_t &A);
    void GetEigenDirectValue(Double_t &ei01, Double_t &ei02, Double_t &ei03);
    void GetEigenValue(Double_t &lambda01, Double_t &lambda02, Double_t &lambda03);
    void GetEigenVector(vector<double> &eigenv01,vector<double> &eigenv02,  vector<double> &eigenv03);
    double ReturnInnerProductOfSquare(TLorentzVector &tlv, vector<double> vec);
    void ReturnTLV(vector<TLorentzVector> &v_tlv, int &v_tlv_size);
    //for check//
    void Check();
    void CheckForEigenVector();

    //////functions used in this program////////////////////
    void CalculateEigenValue();
    void CalculateSafeEigenValue();
    void PrePareSafeMatrix();
    void PrePareMatrix();
    double ReturnXYZ(const int i , const TLorentzVector &tlv);
    void InputEigenValues();
    void InputEigenVectors();
    void SortEigenValues();
    double ReturnPzBoost(TLorentzVector &tlv);
    void GetBoost(TLorentzVector &tlv, double *b);
    void Flip(double &a, double &b){
            double c = a ;
            a=b;
            b=c;
    }

private:
    Int_t m_size;
    vector<TLorentzVector> m_v_tlv;
    int m_v_tlv_size;
    double m_precise;
    TMatrixDSym m_M;
    TVectorD m_eigen;
    TMatrixD m_eigenv;
    double m_lambda01;
    double m_lambda02;
    double m_lambda03;
    vector<double> m_eigenv01;
    vector<double> m_eigenv02;
    vector<double> m_eigenv03;
};
#endif
