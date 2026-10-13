// SPDX-License-Identifier: GPL-3.0-or-later
// Wang et al. (2024), equations 14, 17, 18, with a bounded tie at zero central slope.
#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include "BrineDual.hpp"
template<class Real> inline std::array<Real,3> phaseVector(Real e,const std::array<double,4>& v,
                                      double dx,double dy,Real Q) {
    using saltchem::scalarValue;using saltchem::hypot;
    const double gx=(v[1]-v[0])/(2*dx),gy=(v[3]-v[2])/(2*dy);
    const Real qx=Q*gx,qy=Q*gy;
    Real d[2]={};
    for(int a=0;a<2;a++) {
        Real q=a?qy:qx;double h=a?dy:dx;
        if(scalarValue(q)>0) {
            Real delta=v[2*a+1]-e;
            d[a]=scalarValue(Q*delta)>0?delta/h:Real(0);
        }
        else if(scalarValue(q)<0) {
            Real delta=v[2*a]-e;
            d[a]=scalarValue(Q*delta)>0?-delta/h:Real(0);
        }
        else if(scalarValue(Q)!=0) {
            Real l=(v[2*a]-e)/h,r=(v[2*a+1]-e)/h;
            if(scalarValue(Q)>0){Real high=scalarValue(l)>scalarValue(r)?l:r;d[a]=scalarValue(high)>0?high:Real(0);}
            else {Real low=scalarValue(l)<scalarValue(r)?l:r;d[a]=scalarValue(low)<0?low:Real(0);}
        }
    }
    return {{hypot(d[0],d[1]),qx,qy}};
}
