// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include <deque>
#include <vector>

namespace saltflow {
// Limited-memory type-II Anderson acceleration. MPI reductions make the small
// least-squares problem identical on every rank. It changes iterates only;
// acceptance still uses the original finite-volume equations and mass ledger.
class CouplingAcceleration {
    using Field=Foam::scalarField;
    int depth;
    double beta,lastNorm=0;
    int stalled=0;
    Field lastX,lastF;
    std::deque<Field> dx,df;
    static double dot(const Field& a,const Field& b){double v=Foam::sum(a*b);Foam::reduce(v,Foam::sumOp<double>());return v;}
public:
    CouplingAcceleration(int memory=5,double damping=.8):depth(memory),beta(damping){
        if(depth<0||depth>12||beta<=0||beta>1)throw std::invalid_argument("Invalid Anderson memory/damping");
    }
    Field next(const Field& x,const Field& fixedPoint){
        Field f(fixedPoint-x);const double norm=dot(f,f);
        stalled=lastX.size()&&norm>.9*lastNorm?stalled+1:0;
        const bool restart=stalled>=8;if(restart)stalled=0;
        if(lastX.size()&&!restart&&norm<16*std::max(lastNorm,1e-60)){
            dx.emplace_back(x-lastX);df.emplace_back(f-lastF);
            if(int(dx.size())>depth){dx.pop_front();df.pop_front();}
        }else{dx.clear();df.clear();}
        lastX=x;lastF=f;lastNorm=norm;
        Field step(beta*f);const int n=int(df.size());
        if(n){
            std::vector<std::vector<double>> a(n,std::vector<double>(n+1));double trace=0;
            // All ranks have the same history length. Reduce the small Gram
            // matrix in one message instead of one collective per dot product.
            std::vector<double> products;products.reserve(n*(n+3)/2);
            for(int i=0;i<n;i++){
                for(int j=0;j<=i;j++)products.push_back(Foam::sum(df[i]*df[j]));
                products.push_back(Foam::sum(df[i]*f));
            }
            Foam::reduce(products.data(),int(products.size()),Foam::sumOp<double>(),Foam::UPstream::msgType());
            int k=0;
            for(int i=0;i<n;i++){
                for(int j=0;j<=i;j++)a[i][j]=a[j][i]=products[k++];
                a[i][n]=products[k++];trace+=a[i][i];
            }
            for(int i=0;i<n;i++)a[i][i]+=std::max(trace*1e-12,1e-60);
            bool ok=true;
            for(int i=0;i<n;i++){
                int pivot=i;for(int k=i+1;k<n;k++)if(std::abs(a[k][i])>std::abs(a[pivot][i]))pivot=k;
                std::swap(a[i],a[pivot]);if(std::abs(a[i][i])<1e-60){ok=false;break;}
                for(int k=i+1;k<n;k++){const double r=a[k][i]/a[i][i];for(int j=i;j<=n;j++)a[k][j]-=r*a[i][j];}
            }
            std::vector<double> gamma(n);
            if(ok)for(int i=n-1;i>=0;i--){double b=a[i][n];for(int j=i+1;j<n;j++)b-=a[i][j]*gamma[j];gamma[i]=b/a[i][i];if(!std::isfinite(gamma[i])||std::abs(gamma[i])>1e4)ok=false;}
            if(ok){
                Field proposal(step);for(int i=0;i<n;i++)proposal-=gamma[i]*(dx[i]+beta*df[i]);
                if(dot(proposal,proposal)<=100*std::max(norm,1e-60))step=proposal;
                else{dx.clear();df.clear();}
            }else{dx.clear();df.clear();}
        }
        return x+step;
    }
};
}
