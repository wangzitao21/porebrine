// SPDX-License-Identifier: GPL-3.0-or-later
#include "MineralInterface.hpp"
#include <iostream>
#include <iomanip>
using namespace saltgeom;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double tol,const char* message){require(std::abs(a-b)<=tol*std::max({std::abs(a),std::abs(b),1e-100}),message);}
// Independent steady FV boundary-value solve: inlet concentration 1, retained
// layer(s), several water cells, and a Robin reaction at the receded right end.
double slab(int n,double L,int layers,const ProductLayer& p){
    const double h=L/n,k=.7;
    std::vector<double> base(n),extra(n),g(n+1),lo(n),diag(n),hi(n),rhs(n),c(n);
    for(int i=0;i<n;i++){
        const double x=(i+.5)*h;
        const bool solid=layers==1?(x>=.125&&x<.375):((x>=.125&&x<.25)||(x>=.5&&x<.625));
        const double product=solid?.2:0,water=1-product;
        base[i]=water*water;extra[i]=retainedLayerResistance(water,product,p);
    }
    g[0]=faceDiffusion(base[0],extra[0],extra[0])/(h*.5);
    g[n]=1/(h*.5/faceDiffusion(base.back(),extra.back(),extra.back())+1/k);
    for(int i=1;i<n;i++)g[i]=faceDiffusion(.5*(base[i-1]+base[i]),extra[i-1],extra[i])/h;
    for(int i=0;i<n;i++){diag[i]=g[i]+g[i+1];lo[i]=-g[i];hi[i]=-g[i+1];}
    rhs[0]=g[0];
    for(int i=1;i<n;i++){const double m=lo[i]/diag[i-1];diag[i]-=m*hi[i-1];rhs[i]-=m*rhs[i-1];}
    c.back()=rhs.back()/diag.back();for(int i=n-2;i>=0;i--)c[i]=(rhs[i]-hi[i]*c[i+1])/diag[i];
    const double flux=g[0]*(1-c[0]);
    for(int i=1;i<n;i++)close(g[i]*(c[i-1]-c[i]),flux,2e-10,"Internal face flux is not conserved");
    close(g[n]*c.back(),flux,2e-10,"Surface consumption differs from inlet flux");
    return flux;
}
// A two-dimensional layer with a resolved opening. The opening must provide
// a lateral bypass, rather than applying a global reduction to every sink.
double apertureFlux(bool opening,bool excess,const ProductLayer& p){
    constexpr int nx=32,ny=16;const double h=1./nx;
    std::vector<double> b(nx*ny),r(nx*ny),c(nx*ny,.5),gx((nx+1)*ny),gy(nx*(ny+1));
    for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
        const int k=y*nx+x;const bool gap=opening&&y>=ny/4&&y<3*ny/4;
        const double product=x>=8&&x<12&&!gap?.2:0;
        b[k]=(1-product)*(1-product);r[k]=excess?retainedLayerResistance(1-product,product,p):0;
    }
    for(int y=0;y<ny;y++)for(int x=0;x<=nx;x++){
        const int a=y*nx+std::max(0,x-1),bcell=y*nx+std::min(nx-1,x);
        gx[y*(nx+1)+x]=faceDiffusion(.5*(b[a]+b[bcell]),r[a],r[bcell])/(x==0||x==nx?.5:1);
    }
    for(int y=1;y<ny;y++)for(int x=0;x<nx;x++)gy[y*nx+x]=.5*(b[(y-1)*nx+x]+b[y*nx+x]);
    double change=1;
    for(int it=0;it<20000&&change>1e-13;it++){
        change=0;
        for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
            const int k=y*nx+x;const double w=gx[y*(nx+1)+x],e=gx[y*(nx+1)+x+1],s=gy[y*nx+x],n=gy[(y+1)*nx+x];
            const double next=(w*(x?c[k-1]:1)+e*(x+1<nx?c[k+1]:0)+s*(y?c[k-nx]:0)+n*(y+1<ny?c[k+nx]:0))/(w+e+s+n);
            change=std::max(change,std::abs(next-c[k]));c[k]=next;
        }
    }
    require(change<=1e-13,"2D aperture solve did not converge");
    double inlet=0,outlet=0;
    for(int y=0;y<ny;y++){inlet+=gx[y*(nx+1)]*(1-c[y*nx]);outlet+=gx[y*(nx+1)+nx]*c[y*nx+nx-1];}
    close(inlet,outlet,2e-9,"2D aperture mass balance");
    return inlet/(ny*h);
}
int main(){try{
    std::cout<<std::setprecision(16)<<"test,value\n";
    ProductLayer layer;layer.enabled=layer.crossGrid=true;layer.diffusivity=1;layer.kinetic[2]=.2;
    close(faceDiffusion(2,0,0),2,1e-15,"Zero product changed legacy transport");
    close(faceDiffusion(2,3,8,.25),faceDiffusion(2,8,3,.75),1e-15,"Face orientation changed flux coefficient");
    require(faceDiffusion(2,3,8)>0&&faceDiffusion(2,3,8)<2,"Resistance failed to reduce conductance");
    require(retainedLayerResistance(.8,0,layer)==0,"No-product resistance");
    SurfaceState state;auto& pop=state.population[2];pop.fraction=.02;pop.count=10;pop.supportArea=1e-10;
    saltchem::MineralVector bare{};bare[2]=pop.supportArea;const double V=1e-15;
    pop.retreat=.4*radius(pop,V);
    auto a=measure(state,bare,V,true,layer);
    close(a.attachedFraction[0]+a.residueFraction[0],pop.fraction,1e-14,"Attached and retained product do not partition mass");
    close(a.layerThickness[2]*bare[2]*a.layerCoverage[2]*(1-layer.porosity)/V,a.attachedFraction[0],1e-14,"Retained product counted in local film");
    a=measure(state,saltchem::MineralVector{},V,true,layer);
    require(a.attachedFraction[0]==0&&a.layerThickness[2]==0&&a.capArea[0]>0,"Parent-free product lost independent chemistry");
    close(a.residueFraction[0],pop.fraction,1e-14,"Parent disappearance deleted retained solid");
    require(retainedLayerResistance(.98,a.residueFraction[0],layer)>0,"Parent-free product lost transport resistance");
    ProductLayer sealed=layer;sealed.porosity=0;
    require(faceDiffusion(1,retainedLayerResistance(.8,.2,sealed),0)<1e-95,"Impermeable layer limit");
    // Analytic continuum resistance: clear water + a porous layer occupying
    // product/(1-phi) of the zone + reaction resistance. No face formula used.
    const double layerFraction=.2/(1-layer.porosity);
    const double series=1-layerFraction+layerFraction*layer.tortuosity/layer.porosity;
    const double exact=1/(.75+.25*series+1/.7);
    double previous=1;
    for(int n:{16,32,64,128}){
        const double flux=slab(n,1,1,layer),error=std::abs(flux/exact-1);
        require(error<previous*.55,"Layer-to-receded-interface flux did not converge");previous=error;
        std::cout<<"slab_relative_error_"<<n<<','<<error<<'\n';
    }
    require(previous<.001,"Finest slab flux error");
    const double shortGap=slab(128,1,1,layer),longGap=slab(192,1.5,1,layer);
    close(1/longGap-1/shortGap,.5,2e-11,"Additional intervening water cells did not add gap resistance");
    require(longGap<shortGap,"A longer water gap did not reduce delivery to the receded interface");
    const double two=slab(256,1,2,layer);
    require(std::abs(two/exact-1)<.001,"Two spatially separated layers failed the series-resistance solution");
    const double closed=apertureFlux(false,true,layer),open=apertureFlux(true,true,layer),unblocked=apertureFlux(false,false,layer);
    require(open>closed&&open<1&&unblocked>closed,"Resolved opening failed to provide a 2D bypass");
    std::cout<<"2d_closed_layer_flux,"<<closed<<"\n2d_open_layer_flux,"<<open<<"\n2d_no_excess_flux,"<<unblocked<<'\n';
    std::cout<<"extra_water_gap_resistance,"<<1/longGap-1/shortGap<<"\ntwo_layer_relative_error,"<<std::abs(two/exact-1)<<"\nall_passed,1\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
