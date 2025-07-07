// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include </opt/homebrew/Cellar/boost/1.81.0_1/include/boost/random.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <armadillo>
#include</opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/unsupported/Eigen/CXX11/Tensor>
#define _USE_MATH_DEFINES
#include <chrono>


using namespace Eigen;
using namespace std;
using namespace arma;

int fattoriale(int a){
    int f=1;
    if(a==0){
        f=1;
    }
    if(a>0){
        for(int i=0;i<a;i++){
            f=f*(a-i);
        }
    }
    return f;
}


tuple <mat,int,int> vector_to_matrix(const std::vector< std::vector<int> > &student_vector)
{
  mat y(student_vector.size(),student_vector[0].size());
  
  for(int i = 0; i < student_vector.size(); i++)
  {
    for(int j = 0; j < student_vector[i].size(); j++)
    {
      y(i,j) = student_vector[i][j];
    }
  }
  
  int _n=student_vector.size();
  int _G=student_vector[0].size();
  
  return {y,_n,_G};
}


// [[Rcpp::export]]
int main2(vector<vector<int>> yy){
  mat y;
  int _n;
  int _G;
  tie(y,_n,_G)=vector_to_matrix(yy);
  
    auto start = chrono::steady_clock::now();
    std::vector< std::vector<int> > student_vector;


    int _L,_R,_Q,_iter,_iter0,_iter3,_starting,_burn,_thinn;
    _L=3;
    _R=2;
    _iter=100;
    _iter2=100;
    _starting=6;
    double eps=0.000000001;
    double eps2=0.000000001;
    
    
    _Q = y(0, 1);
    for (int j = 0; j < _n; j++) {
        if (y(j, 1) > _Q) {
            _Q = y(j, 1);
        }
    }

//_vettoreqk contains the n observations ordered according to Q
    vector<int> y0(_n);
    for (int j = 0; j < _n; j++) {
        y0[j]=j+1;
    }
    vector<int> _vettoreqk(_n);
    int xx = 0;
    for (int q = 0; q < _Q; q++) {
        for (int j = 0; j < _n; j++) {
            if (y(j, 0) == q + 1) {
                _vettoreqk[xx] = y0[j];
                xx = xx + 1;
            }
        }
    }

    
    int xxx = 0;
    for (int j = 0; j < _n; j++) {
        if (y0[j] == _vettoreqk[xxx]) {
            _indvecqk[xxx] = j;
            xxx = xxx + 1;
            if (xxx <= _n - 1) {
                j = -1;
            } else if (xxx == _n) {
                j = _n - 1;
            }
        }
    }


    vector<int> _nq(_Q);

    for (int q = 0; q < _Q; q++) {
        _nq[q] = 0;
    }
    for (int q = 0; q < _Q; q++) {
        for (int j = 0; j < _n; j++) {
            if (y(j, 0) == q + 1) {
                _nq[q] += 1;
            }
        }
    }


    int _I = _G-1;
    vector<int> vC(_I);

    vector<int> vC(_I);
    int ma;
    for(int i=0;i<_I;i++){
        ma = y(0,i+2);
        for (int j = 0; j < _n; j++) {
            if (y(j,i+2) > ma) {
                ma = y(j,i+2);
            }
        }
        vC[i]=ma+1;
    }

    int _C; //max between the values that the variables assume +1
    int maa;
    maa = vC[0];
    for (int i = 0; i < _I; i++) {
        if (vC[i] > maa) {
            maa = vC[i];
        }
    }
    _C=maa;

    mat _pi(_I,_L);
    mat _xi(_I, _L);
    mat _w(_n,_L);
    arma::cube _eta(_I,_n,_L);
    arma::cube f(_I,_n,_L);
    vector<double> _pl(_L);
    vector<double> _de(_n);
    vector<double> _de2(_L);
    mat _tau(_n, _L);

    arma::cube _PI(_starting,_I,_L);
    arma::cube _XI(_starting,_I,_L);
    mat _PL(_starting,_L);
    vector<double> LOGL(_starting);

    //------------------------------------------------------------------------
    for (int st = 0; st < _starting; st++) {

        for (int l = 0; l < _L; l++) {
            _pl[l] = (double) 1 / _L;
        }

        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                _pi(i, l) =(double) 0.5;
            }
        }

        mat _xi2(_I, _L,fill::randu);
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                _xi(i, l) =(double) _xi2(i,l);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int i = 0; i < _I; i++) {
                    f(i, j, l) = (double) _pi(i, l) *
                                 (fattoriale(vC[i] - 1) / (fattoriale(y(j, i+1) - 1) * fattoriale(vC[i] - y(j, i+1)))) *
                                 pow((1 - _xi(i, l)), (y(j, i+1) - 1)) * pow(_xi(i, l), (vC[i] - y(j, i+1))) +
                                 (1 - _pi(i, l)) / vC[i];
                }
            }
        }

        mat _f1(_n, _L);
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _f1(j, l) = 1;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int i = 0; i < _I; i++) {
                    _f1(j, l) = (double) _f1(j, l) * f(i, j, l);
                }
            }
        }

        double L, L1;
        vector<double> _flog(_n);
        for (int j = 0; j < _n; j++) {
            _flog[j] = 0;
        }
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _flog[j] += (double) _pl[l] * _f1(j, l);
            }
        }
        for (int j = 0; j < _n; j++) {
            _flog[j] = (double) log(_flog[j]);
        }
        L1 = 0;
        for (int j = 0; j < _n; j++) {
            L1 += (double) _flog[j];
        }

        for (int u = 0; u < _iter; u++) {

            mat num0(_n, _L);
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = (double) log(_f1(j, l)) + log(_pl[l]);
                }
            }
            vector<double> _max0(_n);
            for (int j = 0; j < _n; j++) {
                _max0[j] = num0(j, 0);
                for (int l = 0; l < _L; l++) {
                    if (num0(j, l) > _max0[j]) {
                        _max0[j] = num0(j, l);
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = (double) num0(j, l) - _max0[j];
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    num0(j, l) = (double) exp(num0(j, l));
                }
            }

            vector<double> _den0(_n);
            for (int j = 0; j < _n; j++) {
                _den0[j] = 0;
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _den0[j] += (double) num0(j, l);
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _tau(j, l) = (double) num0(j, l) / _den0[j];
                }
            }

            for (int l = 0; l < _L; l++) {
                _pl[l] = 0;
            }

            for (int l = 0; l < _L; l++) {
                for (int j = 0; j < _n; j++) {
                    _pl[l] += (double) _tau(j, l);
                }
            }

            for (int l = 0; l < _L; l++) {
                _de2[l] = (double) _pl[l];
            }

            for (int l = 0; l < _L; l++) {
                _pl[l] = (double) _pl[l] / _n;
            }


            for (int i = 0; i < _I; i++) {
                for (int j = 0; j < _n; j++) {
                    for (int l = 0; l < _L; l++) {
                        _eta(i, j, l) = (double) (f(i, j, l) - (1 - _pi(i, l)) / vC[i]) / f(i, j, l);
                    }
                }
            }

            mat num1(_I, _L);
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    num1(i, l) = 0;
                }
            }
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int j = 0; j < _n; j++) {
                        num1(i, l) += (double) _tau(j, l) * _eta(i, j, l);
                    }
                }
            }


            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    _pi(i, l) = (double) num1(i, l) / _de2[l];
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                        if (_pi(i, l) >= 0.999) {
                            _pi(i, l) = 0.99;
                        }
                        if (_pi(i, l) < 0.0001) {
                            _pi(i, l) = 0.001;
                        }
                    }
            }

           
            mat me(_I, _L);
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    me(i, l) = 0;
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int j = 0; j < _n; j++) {
                        me(i, l) += (double) _tau(j, l) * _eta(i, j, l) * y(j, i+1);
                    }
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    me(i, l) = (double) me(i, l) / num1(i, l);
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    if(num1(i,l)==0) {
                        me(i, l) = 0;
                    }
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    _xi(i, l) = (double) (vC[i] - me(i, l)) / (vC[i] - 1);
                }
            }


                for (int i = 0; i < _I; i++) {
                    for (int l = 0; l < _L; l++) {
                        if (_xi(i,l) >= 0.999) {
                            _xi(i, l) = 0.99;
                        }
                        if (_xi(i,l) < 0.0001) {
                            _xi(i, l) = 0.001;
                        }
                    }
                }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int i = 0; i < _I; i++) {
                        f(i, j, l) = (double) _pi(i, l) *
                                     (fattoriale(vC[i] - 1) / (fattoriale(y(j, i+1) - 1) * fattoriale(vC[i] - y(j, i+1)))) *
                                     pow((1 - _xi(i, l)), (y(j, i+1) - 1)) * pow(_xi(i, l), (vC[i] - y(j, i+1))) +
                                     (1 - _pi(i, l)) / vC[i];
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _f1(j, l) = 1;
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int i = 0; i < _I; i++) {
                        _f1(j, l) = (double) _f1(j, l) * f(i, j, l);
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                _flog[j] = 0;
            }
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _flog[j] += (double) _pl[l] * _f1(j, l);
                }
            }
            for (int j = 0; j < _n; j++) {
                _flog[j] = (double) log(_flog[j]);
            }
            L = 0;
            for (int j = 0; j < _n; j++) {
                L += (double) _flog[j];
            }
            if(u==_iter-1){
                cout<<"Maximum iteration achieved in the initialization!"<<endl;
            }
            if (abs(L1 - L) < eps) {
                u = _iter - 1;
            }
            L1 = L;

        }

        double min1, temp1;
        double tempp;
        for (int l = 0; l < _L - 1; l++) {
            min1 = l;
            for (int j = l + 1; j < _L; j++) {
                if (_pl[j] < _pl[min1]) {
                    min1 = j;
                }
            }
            temp1 = _pl[l];
            _pl[l] = _pl[min1];
            _pl[min1] = temp1;
            for (int i = 0; i < _I; i++) {
                tempp = _pi(i, l);
                _pi(i, l) = _pi(i, min1);
                _pi(i, min1) = tempp;
            }
            for (int i = 0; i < _I; i++) {
                tempp = _xi(i, l);
                _xi(i, l) = _xi(i, min1);
                _xi(i, min1) = tempp;
            }
        }

        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                _PI(st, i, l) = (double) _pi(i, l);
            }
        }

        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                _XI(st, i, l) = (double) _xi(i, l);
            }
        }

        for (int l = 0; l < _L; l++) {
            _PL(st, l) = (double) _pl[l];
        }

        LOGL[st] = (double) L1;

    }

    double _maxlogl;
    int indicestarting=0;
    _maxlogl = LOGL[0];
    for (int st = 0; st < _starting; st++) {
        if (LOGL[st] > _maxlogl) {
            _maxlogl = LOGL[st];
            indicestarting = st;
        }
    }

    /*cout<<"pi"<<endl;
    for (int i = 0; i < _I; i++) {
        for (int l = 0; l < _L; l++) {
            cout << _PI(indicestarting,i, l) << ",";
        }
        cout<<"\n";
    }
    cout<<"\n";
    cout<<"xi"<<endl;
    for (int i = 0; i < _I; i++) {
        for (int l = 0; l < _L; l++) {
            cout<<_XI(indicestarting,i,l)<<",";
        }
        cout<<"\n";
    }
    cout<<"\n";
    cout<<"pl"<<endl;
    for (int l = 0; l < _L; l++) {
        cout<<_PL(indicestarting,l)<<",";
    }
    cout<<"\n";
    cout<<"logL"<<endl;
    cout<<LOGL[indicestarting]<<endl;*/

    int npar1;
    double BIC1,AIC1,AIC31,CAIC1,SABIC1;
    npar1=2*_I;

    npar1 = _L*npar1+_L-1;
    BIC1 = -2 * LOGL[indicestarting] + log(_n) * npar1;
    SABIC1 = -2 * LOGL[indicestarting] + log((_n+2)/24) * npar1;
    CAIC1 = -2 * LOGL[indicestarting] + (log(_n)+1) * npar1;
    AIC1 = -2 * LOGL[indicestarting] + 2 * npar1;
    AIC31 = -2 * LOGL[indicestarting] + 3 * npar1;

    /*cout << "BIC" << endl;
    cout << BIC1 << endl;
    cout << "AIC" << endl;
    cout << AIC1 << endl;
    cout << "AIC3" << endl;
    cout << AIC31 << endl;
    cout << "CAIC" << endl;
    cout << CAIC1 << endl;
    cout << "SABIC" << endl;
    cout << SABIC1 << endl;*/


    arma::cube _pim(_I,_L,_R);
    arma::cube _xim(_I,_L,_R);
    for (int i = 0; i < _I; i++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _pim(i,l,r)=(double) _PI(indicestarting, i, l);
            }
        }
    }
   
    for (int i = 0; i < _I; i++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _xim(i,l,r)=(double) _XI(indicestarting, i, l);
            }
        }
    }


    mat _z1(_Q, _R);
    vector<double> _numr(_R);
    vector<double> _p_r(_R);
    arma::cube _f(_n, _L, _R);
    mat _pi_l(_L,_R);
    for (int r = 0; r < _R; r++) {
        for (int l = 0; l < _L; l++) {
            _pi_l(l, r) = (double) 1 / _L;
        }
    }
    arma::cube _pf(_n, _L, _R);
    mat _fykq(_n, _R);
    mat _ff2(_Q, _R);
    mat _num3(_Q, _R);
    vector<double> _den3(_Q);
    mat _Prob2(_Q,_R);
    mat _denom(_n, _R);
    arma::cube _Prob3(_n, _L, _R);
    arma::cube _Prob4(_n, _L, _R);
    mat _num4(_L, _R);
    vector<double> _den5(_L);
    vector<double> _den4(_R);
    arma::cube _num5(_I, _L, _R);
    Tensor<double,4> _etam(_I,_n,_L,_R);

    _p_r[0]=0.5;
    while(_p_r[0]==0.5) {
        _z1.zeros(_Q, _R);
        double eqprob2;
        eqprob2 = (double) 1 / _R;
        vector<double> vec3(_R);
        std::fill(vec3.begin(), vec3.end(), eqprob2);
        std::random_device rdtest;
        std::mt19937 generator2(rdtest());
        boost::random::discrete_distribution<int> distribution3(vec3.begin(), vec3.end());
        for (int q = 0; q < _Q; q++) {
            int sample = distribution3(generator2);
            _z1(q, sample) = 1;
        }

        std::fill(_numr.begin(), _numr.end(), 0);

        for (int r = 0; r < _R; r++) {
            for (int q = 0; q < _Q; q++) {
                _numr[r] += _z1(q, r);
            }
        }

        //update p_r

        for (int r = 0; r < _R; r++) {
            _p_r[r] = (double) _numr[r] / _Q;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _f(j, l, r) = 1;
            }
        }
    }

    Tensor<double,4> fm(_I,_n,_L,_R);
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int i = 0; i < _I; i++) {
                    fm(i, j, l, r) = (double) _pim(i, l, r) *
                                 (fattoriale(vC[i] - 1) /
                                  (fattoriale(y(j, i + 1) - 1) * fattoriale(vC[i] - y(j, i + 1)))) *
                                 pow((1 - _xim(i, l, r)), (y(j, i + 1) - 1)) *
                                 pow(_xim(i, l, r), (vC[i] - y(j, i + 1))) +
                                 (1 - _pim(i, l, r)) / vC[i];
                }
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int i = 0; i < _I; i++) {
                    _f(j, l, r) =(double) _f(j, l, r) * fm(i,j,l,r);
                }
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pf(j, l, r) = (double) _pi_l(l, r) * _f(j, l, r);
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            _fykq(j, r) = 0;
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _fykq(j, r) +=(double) _pf(j, l, r);
            }
        }
    }


    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _ff2(q, r) = 1;
        }
    }

    int mm = 0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            for (int j = mm; j < _nq[q] + mm; j++) {
                _ff2(q, r) =(double) _ff2(q, r) * _fykq(_indvecqk[j], r);
            }
        }
        mm = mm + _nq[q];
    }


    //log-liklihood
    double _logL = 0;
    double _logL1 = 0;
    vector<double> _tt(_Q);
    for (int q = 0; q < _Q; q++) {
        _tt[q] = 0;
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _tt[q] +=(double) _p_r[r] * _ff2(q, r);
        }
    }

    for (int q = 0; q < _Q; q++) {
        _logL +=(double) log(_tt[q]);
    }

    double _Lik;
    _Lik=1;
    for (int q = 0; q < _Q; q++) {
        _Lik =(double)_Lik*_tt[q];
    }

    cout<<_logL<<endl;
    cout<<_Lik<<endl;

for (int u1 = 0; u1 < _iter2; u1++) {

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 0;
            }
        }

        mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int j = mm; j < _nq[q] + mm; j++) {
                    _ff2(q, r) +=(double) log(_fykq(_indvecqk[j], r));
                }
            }
            mm = mm + _nq[q];
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) =(double) log(_p_r[r]) + _ff2(q, r);
            }
        }

        vector<double> _maxr(_Q);
        for (int q = 0; q < _Q; q++) {
            _maxr[q] = _num3(q, 0);
            for (int r = 0; r < _R; r++) {
                if (_num3(q, r) > _maxr[q]) {
                    _maxr[q] = _num3(q, r);
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = _num3(q, r) - _maxr[q];
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = exp(_num3(q, r));
            }
        }

        for (int q = 0; q < _Q; q++) {
            _den3[q] = 0;
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _den3[q] += _num3(q, r);
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _Prob2(q, r) = (double) _num3(q, r) / _den3[q];
            }
        }

        for (int r = 0; r < _R; r++) {
            _numr[r] = 0;
        }
        for (int r = 0; r < _R; r++) {
            for (int q = 0; q < _Q; q++) {
                _numr[r] += _Prob2(q, r);
            }
        }

        

        for (int r = 0; r < _R; r++) {
            _p_r[r] = (double) _numr[r] / _Q;
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3(j, l, r) = log(_pi_l(l, r)) + log(_f(j, l, r));
                }
            }
        }

        mat _maxl(_n, _R);
        for (int r = 0; r < _R; r++) {
            for (int j = 0; j < _n; j++) {
                _maxl(j, r) = _Prob3(j, 0, r);
                for (int l = 0; l < _L; l++) {
                    if (_Prob3(j, l, r) > _maxl(j, r)) {
                        _maxl(j, r) = _Prob3(j, l, r);
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3(j, l, r) = _Prob3(j, l, r) - _maxl(j, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3(j, l, r) = exp(_Prob3(j, l, r));
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            _denom[j] = 0;
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _denom[j] += _Prob3(j, l, r);
                }
            }
        }

        mat _denomi(_n, _R);
        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _denomi(j, r) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _denomi(j, r) +=(double) _Prob3(j, l, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _Prob3(j, l, r) = (double) _Prob3(j, l, r) / _denomi(j, r);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int j = 0; j < _n; j++) {
                    _Prob4(j, l, r) = _Prob2(y(j, 0) - 1, r) * _Prob3(j, l, r);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _num4(l, r) = 0;
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int j = 0; j < _n; j++) {
                    _num4(l, r) += _Prob4(j, l, r);
                }
            }
        }

        for (int r = 0; r < _R; r++) {
            _den4[r] = 0;
        }

        for (int r = 0; r < _R; r++) {
            //for (int l = 0; l < _L; l++) {
            for (int j = 0; j < _n; j++) {
                _den4[r] += _Prob2(y(j, 0) - 1, r);
            }
            //}
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                if (_den4[r] != 0) {
                    _pi_l(l, r) = (double) _num4(l, r) / _den4[r];
                }
            }
        }

       
        vector<int> _dephr(_I);
        for (int i = 0; i < _I; i++) {
            _dephr[i] = 0;
        }
        //_dephr[0]=2;
        //_dephr[1]=2;
        //_dephr[2]=2;
        //_dephr[3]=1;
        //_dephr[4]=2;
        //_dephr[5]=2;
        //_dephr[6]=2;
        //_dephr[7]=3;
        for (int i = 0; i < _I; i++) {
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _etam(i, j, l,r) = (double) (fm(i, j, l, r) - (1 - _pim(i, l, r)) / vC[i]) / fm(i, j, l, r);
                    }
                }
            }
        }

        mat _num7(_n, _L);
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _num7(j, l) = 0;
            }
        }
        for (int l = 0; l < _L; l++) {
            for (int j = 0; j < _n; j++) {
                for (int r = 0; r < _R; r++) {
                    _num7(j, l) +=(double) _Prob4(j, l, r);
                }
            }
        }

        for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        _num5(i, l, r) = 0;
                    }
                }
        }
        for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = 0; j < _n; j++) {
                            if (_dephr[i] == 1) {
                                _num5(i, l, r) +=(double) _Prob4(j, l, r)*_etam(i,j,l,r);
                            } else if (_dephr[i] == 0) {
                                _num5(i, l, r) += _num7(j, l)*_etam(i,j,l,r);
                            }
                        }
                    }
                }
        }

        for (int l = 0; l < _L; l++) {
            _den5[l] = 0;
        }

        for (int l = 0; l < _L; l++) {
            for (int j = 0; j < _n; j++) {
                for (int r = 0; r < _R; r++) {
                    _den5[l] += _Prob4(j, l, r);
                }
            }
        }

        for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if (_dephr[i] == 1) {
                            if (_num4(l, r) != 0) {
                                _pim(i, l, r) = (double) _num5(i, l, r) / _num4(l, r);
                            }
                        } else if (_dephr[i] == 0) {
                            if (_den5[l] != 0) {
                                _pim(i, l, r) = (double) _num5(i, l, r) / _den5[l];
                            }
                        }
                    }
                }
        }

        for (int i = 0; i < _I; i++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                if (_pim(i, l, r) >= 0.999) {
                    _pim(i, l, r) = 0.99;
                }
                if (_pim(i, l, r) < 0.0001) {
                    _pim(i, l, r) = 0.001;
                }
            }
        }
    }

        arma::cube mem(_I, _L,_R);
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    mem(i, l,r) = 0;
                }
            }
        }

        for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = 0; j < _n; j++) {
                            if (_dephr[i] == 1) {
                                mem(i, l, r) +=(double) _Prob4(j, l, r)*_etam(i,j,l,r)*y(j,i+1);
                            } else if (_dephr[i] == 0) {
                                mem(i, l, r) +=(double) _num7(j, l)*_etam(i,j,l,r)*y(j,i+1);
                            }
                        }
                    }
                }
        }

        arma::cube _den6(_I,_L,_R);
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _den6(i,l,r)=0;
                }
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    for (int j = 0; j < _n; j++) {
                        _den6(i,l,r)+=(double) _Prob4(j, l, r)*_etam(i,j,l,r);
                    }
                }
            }
        }

        mat _den7(_I,_L);
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                _den7(i, l) = 0;
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    for (int j = 0; j < _n; j++) {
                        _den7(i,l)+=(double) _Prob4(j, l, r)*_etam(i,j,l,r);
                    }
                }
            }
        }

        for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int r = 0; r < _R; r++) {
                        if (_dephr[i] == 1) {
                            if (_den6(i,l, r) != 0) {
                                mem(i, l, r) = (double) mem(i, l, r) / _den6(i,l, r);
                            }
                        } else if (_dephr[i] == 0) {
                            if (_den7(i,l) != 0) {
                                mem(i, l, r) = (double) mem(i, l, r) / _den7(i,l);
                            }
                        }
                    }
                }
        }

        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _xim(i, l, r) = (double) (vC[i] - mem(i, l, r)) / (vC[i] - 1);
                }
            }
        }

        for (int r = 0; r < _R; r++) {
            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    if (_xim(i,l, r) >= 0.999) {
                        _xim(i, l, r) = 0.99;
                    }
                    if (_xim(i,l, r) < 0.0001) {
                        _xim(i, l, r) = 0.001;
                    }
                }
            }
        }

        //log-liklihood
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    _f(j, l, r) = 1;
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    for (int i = 0; i < _I; i++) {
                        fm(i, j, l, r) = (double) _pim(i, l, r) *
                                         (fattoriale(vC[i] - 1) /
                                          (fattoriale(y(j, i + 1) - 1) * fattoriale(vC[i] - y(j, i + 1)))) *
                                         pow((1 - _xim(i, l, r)), (y(j, i + 1) - 1)) *
                                         pow(_xim(i, l, r), (vC[i] - y(j, i + 1))) +
                                         (1 - _pim(i, l, r)) / vC[i];
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                    for (int i = 0; i < _I; i++) {
                        _f(j, l, r) =(double) _f(j, l, r) * fm(i,j,l,r);
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _pf(j, l, r) = (double) _pi_l(l, r) * _f(j, l, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _fykq(j, r) = 0;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _fykq(j, r) +=(double) _pf(j, l, r);
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 1;
            }
        }

        mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int j = mm; j < _nq[q] + mm; j++) {
                    _ff2(q, r) =(double) _ff2(q, r) * _fykq(_indvecqk[j], r);
                }
            }
            mm = mm + _nq[q];
        }

        //log-liklihood
        _logL1 = 0;

        for (int q = 0; q < _Q; q++) {
            _tt[q] = 0;
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _tt[q] +=(double) _p_r[r] * _ff2(q, r);
            }
        }

        for (int q = 0; q < _Q; q++) {
            _logL1 +=(double) log(_tt[q]);
        }

        if(u1==_iter2-1){
            cout<<"Maximum iteration achieved in MCUB"<<endl;
        }

        if (abs(_logL - _logL1) < eps2) {
            u1 = _iter2 - 1;
        }
        _logL = _logL1;

    }

    cout<<"pr"<<endl;
    for(int r=0;r<_R;r++){
        cout<<_p_r[r]<<",";
    }
    cout<<"\n";

    cout << "_pi_l(l,0)";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        cout << _pi_l(l, 0) << ",";
    }
    cout << "\n";

    cout << "_pi_l(l,1)";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        cout << _pi_l(l, 1) << ",";
    }
    cout << "\n";

    cout << "_pi(0,l,r):";
    cout << "\n";
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
            cout << _pim(0, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_pi(1,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _pim(1, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_pi(2,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _pim(2, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_pi(3,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _pim(3, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_pi(4,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _pim(4, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_xi(0,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _xim(0, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_xi(1,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _xim(1, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_xi(2,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _xim(2, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_xi(3,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _xim(3, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    cout << "_xi(4,l,r):";
    cout << "\n";
    for (int l = 0; l < _L; l++) {
        for (int r = 0; r < _R; r++) {
            cout << _xim(4, l, r) << ",";
        }
        cout << '\n';
    }
    cout << "\n";

    /*double CL;
    CL=0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            CL+=_Prob2(q,r)*log(_p_r[r]);
        }
    }
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                CL += _Prob4(j,l, r) * log(_pi_l(l,r)*_f(j,l,r));
            }
        }
    }

    cout<<CL<<endl;

    double EN,EN1,EN2;
    EN1=0;
    EN2=0;
    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            EN1+=-_Prob2(q,r)*log(_Prob2(q,r));
        }
    }
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                EN2 += -_Prob4(j,l, r) * log(_Prob3(j,l,r));
            }
        }
    }
    cout<<EN1+EN2<<endl;
    cout<<_logL1<<endl;
    cout<<CL+EN1+EN2<<endl;

    int npar;
    npar=_L*2*_I;
    double AIC3e,AICe,BICe,BICge,SABICe,SABICge;
    AIC3e=-2*CL+ 3 * (npar+_R-1+(_R*(_L-1)));
    AICe=-2*CL+ 2 * (npar+_R-1+(_R*(_L-1)));
    BICe = -2 * CL + log(_n) * (npar+_R-1+(_R*(_L-1)));
    BICge = -2 * CL + log(_Q) * (npar+_R-1+(_R*(_L-1)));
    SABICe = -2 * CL + log((_n+2)/24) * (npar+_R-1+(_R*(_L-1)));
    SABICge = -2 * CL + log((_Q+2)/24) * (npar+_R-1+(_R*(_L-1)));

    cout<<"AIC3e"<<endl;
    cout<<AIC3e<<endl;
    cout<<"AICe"<<endl;
    cout<<AICe<<endl;
    cout<<"BICe"<<endl;
    cout<<BICe<<endl;
    cout<<"BICge"<<endl;
    cout<<BICge<<endl;
    cout<<"SABICe"<<endl;
    cout<<SABICe<<endl;
    cout<<"SABICge"<<endl;
    cout<<SABICge<<endl;

    double AIC3,AIC,BIC,BICg,SABIC,SABICg;
    AIC3=-2*(CL+EN1+EN2)+ 3 * (npar+_R-1+(_R*(_L-1)));
    AIC=-2*(CL+EN1+EN2)+ 2 * (npar+_R-1+(_R*(_L-1)));
    BIC = -2 * (CL+EN1+EN2) + log(_n) * (npar+_R-1+(_R*(_L-1)));
    BICg = -2 * (CL+EN1+EN2) + log(_Q) * (npar+_R-1+(_R*(_L-1)));
    SABIC = -2 * (CL+EN1+EN2) + log((_n+2)/24) * (npar+_R-1+(_R*(_L-1)));
    SABICg = -2 * (CL+EN1+EN2) + log((_Q+2)/24) * (npar+_R-1+(_R*(_L-1)));

    cout<<"AIC3"<<endl;
    cout<<AIC3<<endl;
    cout<<"AIC"<<endl;
    cout<<AIC<<endl;
    cout<<"BIC"<<endl;
    cout<<BIC<<endl;
    cout<<"BICg"<<endl;
    cout<<BICg<<endl;
    cout<<"SABIC"<<endl;
    cout<<SABIC<<endl;
    cout<<"SABICg"<<endl;
    cout<<SABICg<<endl;*/


    //classification
    mat _Prob2e(_n,_R);
    for (int j = 0; j < _n; j++) {
        for (int r = 0; r < _R; r++) {
            _Prob2e(j, r) = _Prob2(y(j, 1) - 1, r);
        }
    }
    vector<int> _indr(_n);
    std::fill(_indr.begin(),_indr.end(),0);
    for(int j=0;j<_n;j++){
        double _max2=_Prob2e(j,0);
        for(int r=0;r<_R;r++){
            if(_Prob2e(j,r)>_max2){
                _max2=_Prob2e(j,r);
                _indr[j]=r;
            }
        }
    }
    vector<int> _indL(_n);
    for (int j = 0; j < _n; j++) {
        double _maxl = _Prob3(j, 0,_indr[j]);
        for (int l = 0; l < _L; l++) {
            if (_Prob3(j, l,_indr[j]) > _maxl) {
                _maxl = _Prob3(j, l,_indr[j]);
                _indL[j]=l;
            }
        }
    }

    vector<int> _indR(_Q);
    std::fill(_indR.begin(),_indR.end(),0);
    for(int q=0;q<_Q;q++){
        double _max2=_Prob2(q,0);
        for(int r=0;r<_R;r++){
            if(_Prob2(q,r)>_max2){
                _max2=_Prob2(q,r);
                _indR[q]=r;
            }
        }
    }

    /*cout<<"cla2"<<endl;
    for (int q = 0; q < _Q; q++) {
        cout<<_indR[q]<<",";
    }
    cout<<"\n";
    cout<<"cla1"<<endl;
    for (int j = 0; j < _n; j++) {
        cout<<_indL[j]<<",";
    }
    cout<<"\n";*/

    auto end = chrono::steady_clock::now();
    int elapsed_time = chrono::duration_cast<chrono::seconds>(end - start).count();
    std::cout << "Time: " << elapsed_time << " sec" << std::endl;

}