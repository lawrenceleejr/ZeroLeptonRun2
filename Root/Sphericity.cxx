#include "ZeroLeptonRun2/Sphericity.h"


void Sphericity::ReturnTLV(vector<TLorentzVector> &v_tlv, int &v_tlv_size){
    v_tlv = m_v_tlv;
    v_tlv_size = m_v_tlv_size;
}
double Sphericity::ReturnInnerProductOfSquare(TLorentzVector &tlv, vector<double> vec){
   return fabs(
        tlv.Px()*vec.at(0)+
        tlv.Py()*vec.at(1)+
        tlv.Pz()*vec.at(2)
   );
}
void Sphericity::GetEigenDirectValue(Double_t &ei01, Double_t &ei02, Double_t &ei03){
    ei01 = m_eigen[0];
    ei02 = m_eigen[1];
    ei03 = m_eigen[2];
}
void Sphericity::GetEigenValue(Double_t &lambda01, Double_t &lambda02, Double_t &lambda03){
    lambda01 = m_lambda01;
    lambda02 = m_lambda02;
    lambda03 = m_lambda03;
}
void Sphericity::GetEigenVector(vector<double> &eigenv01,vector<double> &eigenv02,  vector<double> &eigenv03){
    eigenv01 = m_eigenv01;
    eigenv02 = m_eigenv02;
    eigenv03 = m_eigenv03;
}
void Sphericity::SetTLV(const vector<TLorentzVector> &v_tlv, int v_tlv_size){
    m_v_tlv_size = v_tlv_size;
    m_v_tlv = v_tlv;
}
void Sphericity::GetBoostedSphericity(Double_t &S, double &ST, Double_t &A){
    TLorentzVector tlv;tlv.SetPtEtaPhiM(0.,0.,0.,0.);
    for(int i =0; i < m_v_tlv_size; i++){
       tlv += m_v_tlv.at(i); 
    }
    double  b = ReturnPzBoost(tlv);
    for(int i =0; i < m_v_tlv_size; i++){
       m_v_tlv.at(i).Boost(0.,0.,b);
    }
    CalculateEigenValue();
    S = (m_lambda02+m_lambda03)*3./2.;
    ST = ( 0. == (m_lambda01+m_lambda02))? 0. : 2.*m_lambda02/(m_lambda01+m_lambda02);
    A = 3.*m_lambda03/2.;
}
double Sphericity::ReturnPzBoost(TLorentzVector &tlv){ 
    return (0. != tlv.E())? -tlv.Pz()/tlv.E() : 0. ;
}
void Sphericity::GetSphericity(Double_t &S, Double_t &ST, Double_t &A){
    CalculateEigenValue();
    S = (m_lambda02+m_lambda03)*3./2.;
    ST = ( 0. == (m_lambda01+m_lambda02))? 0. : 2.*m_lambda02/(m_lambda01+m_lambda02);
    A = 3.*m_lambda03/2.;
}
void Sphericity::GetSafeSphericity(Double_t &S, Double_t &ST, Double_t &A){
    CalculateSafeEigenValue();
    S = (m_lambda02+m_lambda03)*3./2.;
    ST = ( 0. == (m_lambda01+m_lambda02))? 0. : 2.*m_lambda02/(m_lambda01+m_lambda02);
    A = 3.*m_lambda03/2.;
}
void Sphericity::CalculateSafeEigenValue(){
    PrePareSafeMatrix();
    TMatrixDSymEigen EM(m_M);
    //get eigen values//
    m_eigen.ResizeTo(m_size);
    m_eigen = EM.GetEigenValues();
    //get eigen vectors//
    m_eigenv.ResizeTo(m_size, m_size);
    m_eigenv = EM.GetEigenVectors();

    InputEigenVectors();
    InputEigenValues();
    SortEigenValues();
}
void Sphericity::CalculateEigenValue(){
    PrePareMatrix();
    TMatrixDSymEigen EM(m_M);
    //get eigen values//
    m_eigen.ResizeTo(m_size);
    m_eigen = EM.GetEigenValues();
    //get eigen vectors//
    m_eigenv.ResizeTo(m_size, m_size);
    m_eigenv = EM.GetEigenVectors();

    InputEigenVectors();
    InputEigenValues();
    SortEigenValues();
}
void Sphericity::PrePareSafeMatrix(){
    m_size =3;
    m_M.ResizeTo(m_size, m_size);
    TArrayD tmp(m_size*m_size);

    for( Int_t i =0; i < m_size*m_size; i++) tmp[i]=0.;   
    m_M = TMatrixDSym(m_size, tmp.GetArray()); 
    TMatrixDSym tmp_M(m_size, tmp.GetArray());
    double p,q;
    for(int k =0; k < m_v_tlv_size; k++){
        for(int i=0; i < 3; i++ ){
            for(int j =i; j < 3; j++){
                    p = ReturnXYZ(i, m_v_tlv.at(k));
                    q = ReturnXYZ(j, m_v_tlv.at(k));
                    tmp[3*i+j]=tmp[i+3*j]=p*q/m_v_tlv.at(k).P();
            }
        }
        tmp_M.SetMatrixArray(tmp.GetArray());
        m_M+=tmp_M;
    }
}
void Sphericity::PrePareMatrix(){
    m_size =3;
    m_M.ResizeTo(m_size, m_size);
    TArrayD tmp(m_size*m_size);

    for( Int_t i =0; i < m_size*m_size; i++) tmp[i]=0.;   
    m_M = TMatrixDSym(m_size, tmp.GetArray()); 
    TMatrixDSym tmp_M(m_size, tmp.GetArray());
    double p,q;
    for(int k =0; k < m_v_tlv_size; k++){
        for(int i=0; i < 3; i++ ){
            for(int j =i; j < 3; j++){

                    p = ReturnXYZ(i, m_v_tlv.at(k));
                    q = ReturnXYZ(j, m_v_tlv.at(k));
                    tmp[3*i+j]=tmp[i+3*j]=p*q;
            }
        }
        tmp_M.SetMatrixArray(tmp.GetArray());
        m_M+=tmp_M;
    }
}
double Sphericity::ReturnXYZ(const int i , const TLorentzVector &tlv){
    switch(i){
        case 0:
            return tlv.Px();
        case 1:
            return tlv.Py();
        case 2:
            return tlv.Pz();
        default:
            return 0.;
    }
}
void Sphericity::InputEigenValues(){
    for(int i =0; i < m_eigen.GetNoElements(); i++){
        switch(i){
            case 0 : 
                m_lambda01 = m_eigen[0];
                break;
            case 1 : 
                m_lambda02 = m_eigen[1];
                break;
            case 2 : 
                m_lambda03 = m_eigen[2];
                break;
            default :
                break;
        } 
    }
    double scale = m_lambda01+m_lambda02+m_lambda03; 
    m_lambda01/=scale;
    m_lambda02/=scale;
    m_lambda03/=scale;    
}
void Sphericity::SortEigenValues(){
    if(m_lambda02 > m_lambda01){
        m_eigenv01.swap(m_eigenv02);
        Flip(m_lambda01, m_lambda02);
    }
    if(m_lambda03 > m_lambda02){
        m_eigenv02.swap(m_eigenv03);
        Flip(m_lambda02, m_lambda03);
    }
    if(m_lambda02 > m_lambda01){
        m_eigenv01.swap(m_eigenv02);
        Flip(m_lambda01, m_lambda02);
    }
}
void Sphericity::Check(){
    PrePareMatrix();
    #if 0
    m_M.Print();
    #endif 
    TMatrixDSymEigen EM(m_M);
    m_eigen.ResizeTo(m_size);
    m_eigen = EM.GetEigenValues();
    //get eigen vectors//
    m_eigenv.ResizeTo(m_size, m_size);
    m_eigenv = EM.GetEigenVectors();
    InputEigenVectors();
    InputEigenValues();
    cout << "eigenvalues are " << m_lambda01 << ", " << m_lambda02 << ", " << m_lambda03 << endl;
    SortEigenValues();
    cout << "sorted eigenvalues are " << m_lambda01 << ", " << m_lambda02 << ", " << m_lambda03 << endl;
    m_eigenv.Print();
    for(int i=0; i < 3 ; i ++){
        cout << m_eigenv01.at(i) << " " << m_eigenv02.at(i) << " " << m_eigenv03.at(i) << endl;
    }
   
    double S = (m_lambda02+m_lambda03)*3./2.;
    double ST = ( 0. == (m_lambda01+m_lambda02))? 0. : 2.*m_lambda02/(m_lambda01+m_lambda02);
    double A = 3.*m_lambda03/2.;
    cout << "S = " << S << "ST = " << ST << "A = " << A <<endl;    
}
void Sphericity::CheckForEigenVector(){
    PrePareMatrix();
    #if 0
    m_M.Print();
    #endif 
    TMatrixDSymEigen EM(m_M);
    m_eigen.ResizeTo(m_size);
    m_eigen = EM.GetEigenValues();
//    InputEigenValues();
    cout << "eigenvalues are " << m_eigen[0] << ", " << m_eigen[1] << ", " << m_eigen[2] << endl;
    m_eigenv.ResizeTo(m_size, m_size);
    m_eigenv = EM.GetEigenVectors();
    m_eigenv.Print();
    TMatrixD M(m_M);
    M*=m_eigenv;
    TMatrixD tmp;
    tmp.ResizeTo(m_size, m_size);
    tmp.TMult(m_eigenv, M);
    tmp.Print();
    const Double_t *array = tmp.GetMatrixArray();
    for(int i =0; i < 3; i++){
        for(int j =0; j < 3; j++){
            cout << " (" << i << "," << j << ") = " << array[i*3+j];
        }
        cout << endl;
    }
//    SortEigenValues();
//    cout << "sorted eigenvalues are " << m_lambda01 << ", " << m_lambda02 << ", " << m_lambda03 << endl;
}
void Sphericity::InputEigenVectors(){
    const Double_t *array = m_eigenv.GetMatrixArray();
    vector<vector<double> > tmp_eigenv;
    tmp_eigenv.resize(3);
    for(int i =0; i < 3; i++){
        tmp_eigenv.at(i).resize(3);
        for(int j =0; j < 3; j++){
	  //cout << " (" << i << "," << j << ") = " << array[i*3+j];
            tmp_eigenv.at(i).at(j) = array[i*3+j];
        }
        //cout << endl;
    }
    for(int i =0; i < 3 ;i++){
        m_eigenv01.push_back(tmp_eigenv.at(i).at(0));
        m_eigenv02.push_back(tmp_eigenv.at(i).at(1));
        m_eigenv03.push_back(tmp_eigenv.at(i).at(2));
    }
//    SortEigenValues();
//    cout << "sorted eigenvalues are " << m_lambda01 << ", " << m_lambda02 << ", " << m_lambda03 << endl;
}

