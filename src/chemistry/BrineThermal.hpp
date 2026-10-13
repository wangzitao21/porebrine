// SPDX-License-Identifier: GPL-3.0-or-later
// PHREEQC 3: k_calc, calc_rho_0, calc_dielectrics (USGS public domain).
#pragma once
#include "BrineSystem.hpp"
#include <algorithm>
namespace saltchem {
struct EquilibriumConstant {
    double reference=0,enthalpyKJ=0;
    std::array<double,6> coefficients{};
    bool analytic=false;
    double at(double T)const {
        if(analytic){const auto& a=coefficients;return a[0]+a[1]*T+a[2]/T+a[3]*std::log10(T)+a[4]/(T*T)+a[5]*T*T;}
        return reference-enthalpyKJ*(298.15-T)/(ln10*.00831470*T*298.15);
    }
};
struct WaterProperties {
    template<class Real> struct StateType {Real rho,aphi,dhb,av,qBorn,pressure;};
    using State=StateType<double>;
    double temperature=298.15,pressure=1,waterT=298.15;
    double rhoSat=0,p0=0,p1=0,p2=0,p3=0,pSat=0;
    double dielectricC=0,dielectricB=0,dielectric1000=0;
    double e2=0,epsilonDerivative=0,qBorn=0,dielectricPressure=0;
    template<class Real> void dielectric(Real pa,Real& e,Real& derivative,Real& q)const {
        const Real pb=pa*1.01325;
        const Real epsilon=dielectric1000+dielectricC*log((dielectricB+pb)/(dielectricB+1000));
        if(!(scalarValue(epsilon)>0)||!std::isfinite(scalarValue(epsilon)))throw std::domain_error("Invalid water dielectric state");
        e=1.671008e-3/(epsilon*waterT);
        derivative=dielectricC/(dielectricB+pb)*1.01325/epsilon;
        q=dielectricC/(dielectricB+pb)/(epsilon*epsilon)*41.84004;
    }
    void set(double T,double pa){
        if(!std::isfinite(T)||T<=228||!std::isfinite(pa)||pa<=0)
            throw std::invalid_argument("Temperature must exceed 228 K and absolute pressure must be positive");
        const bool thermalChanged=T!=temperature||rhoSat==0;
        temperature=T;pressure=pa;waterT=std::min(T,623.15);
        if(thermalChanged){
        const double tc=waterT-273.15,th=1-waterT/647.096;
        rhoSat=322*(1+1.99274064*std::pow(th,1./3)+1.09965342*std::pow(th,2./3)-.510839303*std::pow(th,5./3)-1.75493479*std::pow(th,16./3)-45.5170352*std::pow(th,43./3)-6.7469445e5*std::pow(th,110./3));
        p0=.05188+tc*(-4.1885519e-4+tc*(6.6780748e-6+tc*(-3.6648699e-8+tc*8.3501912e-11)));
        p1=-6.0251348e-6+tc*(3.6696407e-7+tc*(-9.2056269e-9+tc*(6.7024182e-11-tc*1.5947241e-13)));
        p2=-2.2983596e-9+tc*(-4.0133819e-10+tc*(1.2619821e-11+tc*(-9.8952363e-14+tc*2.3363281e-16)));
        p3=7.0517647e-11+tc*(6.8566831e-12+tc*(-2.2829750e-13+tc*(1.8113313e-15-tc*4.2475324e-18)));
        pSat=std::exp(11.6702-3816.44/(waterT-46.13));
        dielectricC=-2.0525+3115.9/(waterT-182.89);
        dielectricB=-8032.5+4214200/waterT+2.1417*waterT;
        dielectric1000=342.79*std::exp(waterT*(-.0050866+waterT*9.469e-7));
        }
        dielectricPressure=std::max(pa,pSat);dielectric(dielectricPressure,e2,epsilonDerivative,qBorn);
    }
    template<class Real> StateType<Real> at(Real activity)const {
        const Real saturation=pSat*(scalarValue(activity)<1?activity:Real(1));
        const Real pa=scalarValue(saturation)>pressure?saturation:Real(pressure),p=pa-saturation+1e-6;
        Real rho=rhoSat+p*(p0+p*(p1+p*(p2+sqrt(p)*p3)));
        if(scalarValue(rho)<.01)rho=Real(.01);
        const Real kappa=(p0+p*(2*p1+p*(3*p2+sqrt(p)*3.5*p3)))/rho;
        rho/=1000;
        Real e=e2,derivative=epsilonDerivative,q=qBorn;
        if(scalarValue(pa)!=dielectricPressure)dielectric(pa,e,derivative,q);
        const Real dh=sqrt(8*3.14159265358979323846*6.02252e23*e*rho/1000);
        return {rho,dh*e/6,dh/1e8,dh*e*.0820597*1000*waterT*(derivative-kappa/3),q,pa};
    }
};
}
