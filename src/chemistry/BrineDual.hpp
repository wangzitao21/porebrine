// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
namespace saltchem {
using std::exp;using std::expm1;using std::log;using std::log10;using std::log1p;
using std::sqrt;using std::pow;using std::tanh;using std::hypot;using std::abs;
template<int N> struct Dual {
    double value=0;std::array<double,N> derivative{};
    Dual()=default;Dual(double v):value(v){}
    static Dual variable(double v,int i){Dual x(v);x.derivative[i]=1;return x;}
    Dual& operator+=(const Dual& b){value+=b.value;for(int i=0;i<N;i++)derivative[i]+=b.derivative[i];return *this;}
    Dual& operator-=(const Dual& b){value-=b.value;for(int i=0;i<N;i++)derivative[i]-=b.derivative[i];return *this;}
    Dual& operator*=(const Dual& b){for(int i=0;i<N;i++)derivative[i]=derivative[i]*b.value+value*b.derivative[i];value*=b.value;return *this;}
    Dual& operator/=(const Dual& b){const double inv=1/b.value;value*=inv;for(int i=0;i<N;i++)derivative[i]=(derivative[i]-value*b.derivative[i])*inv;return *this;}
};
inline double scalarValue(double x){return x;}
template<int N> double scalarValue(const Dual<N>& x){return x.value;}
template<int N> Dual<N> operator+(Dual<N> a,const Dual<N>& b){return a+=b;}
template<int N> Dual<N> operator-(Dual<N> a,const Dual<N>& b){return a-=b;}
template<int N> Dual<N> operator*(Dual<N> a,const Dual<N>& b){return a*=b;}
template<int N> Dual<N> operator/(Dual<N> a,const Dual<N>& b){return a/=b;}
template<int N> Dual<N> operator-(Dual<N> a){a.value=-a.value;for(double& d:a.derivative)d=-d;return a;}
#define SALT_DUAL_SCALAR(OP) \
template<int N> Dual<N> operator OP(Dual<N> a,double b){return a OP Dual<N>(b);} \
template<int N> Dual<N> operator OP(double a,Dual<N> b){return Dual<N>(a) OP b;}
SALT_DUAL_SCALAR(+) SALT_DUAL_SCALAR(-) SALT_DUAL_SCALAR(*) SALT_DUAL_SCALAR(/)
#undef SALT_DUAL_SCALAR
template<int N> Dual<N> chain(const Dual<N>& x,double value,double slope){Dual<N> r(value);for(int i=0;i<N;i++)r.derivative[i]=slope*x.derivative[i];return r;}
template<int N> Dual<N> exp(const Dual<N>& x){const double v=std::exp(x.value);return chain(x,v,v);}
template<int N> Dual<N> expm1(const Dual<N>& x){return chain(x,std::expm1(x.value),std::exp(x.value));}
template<int N> Dual<N> log(const Dual<N>& x){if(!(x.value>0))throw std::domain_error("Derivative log domain");return chain(x,std::log(x.value),1/x.value);}
template<int N> Dual<N> log10(const Dual<N>& x){return log(x)/std::log(10.);}
template<int N> Dual<N> log1p(const Dual<N>& x){return chain(x,std::log1p(x.value),1/(1+x.value));}
template<int N> Dual<N> sqrt(const Dual<N>& x){
    if(x.value==0){for(double d:x.derivative)if(d!=0)throw std::domain_error("Singular square-root derivative");return Dual<N>(0);}
    const double v=std::sqrt(x.value);return chain(x,v,.5/v);
}
template<int N> Dual<N> pow(const Dual<N>& x,double a){
    if(a==0)return Dual<N>(1);
    if(x.value==0){if(a<1)throw std::domain_error("Singular power derivative");return chain(x,0,a==1?1:0);}
    const double v=std::pow(x.value,a);return chain(x,v,a*v/x.value);
}
template<int N> Dual<N> abs(const Dual<N>& x){return x.value<0?-x:x;}
template<int N> Dual<N> tanh(const Dual<N>& x){const double v=std::tanh(x.value);return chain(x,v,1-v*v);}
template<int N> Dual<N> hypot(const Dual<N>& a,const Dual<N>& b){
    const double v=std::hypot(a.value,b.value);Dual<N> r(v);
    if(v>0)for(int i=0;i<N;i++)r.derivative[i]=(a.value*a.derivative[i]+b.value*b.derivative[i])/v;
    return r;
}
}
