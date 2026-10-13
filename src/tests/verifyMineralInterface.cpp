// SPDX-License-Identifier: GPL-3.0-or-later
#include "MineralInterface.hpp"
#include "BrineReactions.hpp"
#include "BrineInlet.hpp"
#include <iostream>
#include <iomanip>
using namespace saltgeom;
using namespace saltchem;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double relative,const char* message){require(std::abs(a-b)<=relative*std::max({std::abs(a),std::abs(b),1e-100}),message);}
std::array<double,8> conserved(const Inventory& s){std::array<double,8> a{};for(int i=0;i<7;i++)a[i]=s.aqueous[i];a[7]=s.waterEquivalent;for(int j=0;j<7;j++){for(int i=0;i<7;i++)a[i]+=stoichiometry[j][i]*s.mineral[j];a[7]+=hydration[j]*waterMolarMass*s.mineral[j];}return a;}
int main(int argc,char** argv){
    try{
        require(argc==2,"Supply pitzer.dat");std::cout<<std::setprecision(16)<<"test,value\n";
        constexpr double ep=1e-5,h=1e-5,v=1e-15;
        std::array<std::array<double,4>,7> neighbours{};std::array<Point,7> normals{};normals[2]={1,0};normals[0]={1,0};
        MineralVector f{};f[2]=.4*(1-ep);
        auto a=reconstruct(f,0,neighbours,normals,h,h,h,ep);
        close(a.solid[2].area(),.4*h*h,1e-12,"PLIC volume");close(a.internalArea[2],h*h,1e-12,"Flat PLIC area");
        f[0]=.1*(1-ep);auto coated=reconstruct(f,0,neighbours,normals,h,h,h,ep);
        require(coated.internalArea[2]==0,"Buried primary surface is reactive");close(coated.internalArea[0],h*h,1e-12,"Outer coating area");
        // Equal-volume mineral replacement changes material interfaces while
        // leaving porosity and total solid volume exactly unchanged.
        f[2]=.3*(1-ep);f[0]=.2*(1-ep);auto replaced=reconstruct(f,0,neighbours,normals,h,h,h,ep);
        close(coated.liquid.area(),replaced.liquid.area(),1e-12,"Equal-volume replacement changed pore space");
        std::cout<<"flat_surface_and_replacement,1\n";
        // Isolated hemispherical crystal: V=2*pi*r^3/3, A=2*pi*r^2.
        SurfaceState single;const double r=1e-7;single.population[2]={2*pi*r*r*r/(3*v),1};MineralVector bare{};bare[2]=h*h;
        auto measured=measure(single,bare,v);close(measured.capArea[0],2*pi*r*r,1e-12,"Single crystal growth area");
        close(measured.coverage[2],pi*r*r/(h*h),1e-12,"Single footprint area");
        single.population[2].fraction*=8;measured=measure(single,bare,v);close(measured.capArea[0],8*pi*r*r,1e-12,"Area-volume scaling");
        auto merged=single;merged.population[2].fraction*=1e4;merged.population[2].count=1e4;
        const double stored=subgrid(merged)[0];measured=measure(merged,bare,v);
        require(measured.coverage[2]>0&&measured.coverage[2]<=1&&measured.capArea[0]<=2*bare[2],"Impingement bounds");close(subgrid(merged)[0],stored,1e-15,"Impingement lost volume");
        std::cout<<"single_crystal_area_and_impingement,1\n";
        ProductLayer layer;layer.enabled=true;layer.porosity=.5;layer.kinetic[2]=.2;layer.validate();
        close(layer.transmission(1e-6,2),.625,1e-14,"Film resistance-in-series solution");
        close(layer.transmission(0,2),1,1e-15,"Vanishing film resistance");
        auto impermeable=layer;impermeable.porosity=0;
        require(impermeable.transmission(1e-6,2)==0,"Impermeable film limit");
        auto fastDiffusion=layer;fastDiffusion.diffusivity=1e6;
        close(fastDiffusion.transmission(1e-6,2),1,1e-14,"Infinite diffusion limit");
        SurfaceState porous;auto& patch=porous.population[2];
        patch.fraction=.1;patch.count=400;patch.supportArea=bare[2];
        const auto fresh=measure(porous,bare,v,true,layer);
        require(fresh.coverage[2]>.9&&fresh.accessibility[2]>.1&&fresh.accessibility[2]<1,"Porous coating became a zero-flux mask");
        auto tinySupport=bare;tinySupport[2]*=.001;
        const auto tiny=measure(porous,tinySupport,v,true,layer);
        close(tiny.coverage[2],fresh.coverage[2],1e-14,"Shrinking substrate amplified contact coverage");
        close(tiny.accessibility[2],fresh.accessibility[2],1e-14,"Shrinking substrate artificially thickened the film");
        close(tiny.capArea[0],fresh.capArea[0],1e-14,"Independent product area shrank with parent");
        MineralVector parentBefore{},parentAfter{},cap= fresh.capArea;
        parentBefore[0]=parentAfter[0]=patch.fraction;parentBefore[2]=.4;parentAfter[2]=.3;
        auto parentMeasure=fresh;parentMeasure.exposed[2]=bare[2];
        grow(porous,parentBefore,parentAfter,parentMeasure,cap,v,3.5e-6,normals,bare,layer,ep);
        auto separated=measure(porous,bare,v,true,layer);
        close(separated.coverage[2],fresh.coverage[2],1e-14,"Normal recession bypassed an unresolved local film");
        close(separated.attachedFraction[0],.1,1e-15,"Normal recession prematurely transferred a film to spatial residue");
        close(separated.capArea[0],fresh.capArea[0],1e-14,"Recession changed independent product area");
        require(separated.layerCoverage[2]>.9&&separated.accessibility[2]>0,"Receded local layer lost its finite diffusion resistance");
        parentBefore=parentAfter;parentAfter[2]=0;
        parentMeasure=separated;parentMeasure.exposed[2]=bare[2];
        grow(porous,parentBefore,parentAfter,parentMeasure,cap,v,3.5e-6,normals,bare,layer,ep);
        MineralVector absentParent{};updatePartition(porous,absentParent,normals);const auto retained=measure(porous,absentParent,v,true,layer);
        close(subgrid(porous)[0],.1,1e-15,"Parent depletion changed product inventory");
        require(retained.exposed[0]>0&&retained.coverage[2]==0,"Parent depletion disabled independent product reactions");
        auto halfDissolved=parentAfter;halfDissolved[0]*=.5;
        grow(porous,parentAfter,halfDissolved,retained,retained.capArea,v,3.5e-6,normals,absentParent,layer,ep);
        close(subgrid(porous)[0],.05,1e-15,"Retained product dissolution did not reduce its population volume");
        std::cout<<"porous_film_limits_and_reference_support,1\nrecession_partition_and_independent_product_reaction,1\n";
        auto promoted=single;MineralVector before{},after{},total{};before[0]=after[0]=subgrid(single)[0];before[2]=after[2]=.4;total[0]=measure(single,bare,v).capArea[0];
        grow(promoted,before,after,measure(single,bare,v),total,v,r,normals);
        require(subgrid(promoted)[0]==0&&promoted.promotions[0]==1,"Resolved transition failed");close(after[0],before[0],1e-15,"Promotion changed mass");
        std::cout<<"conservative_promotion,1\n";
        auto exactSubgrid=single;MineralVector grown=after;
        for(int step=0;step<20;step++){
            auto prior=grown;grown[0]*=1.17;
            const auto current=measure(exactSubgrid,bare,v);
            grow(exactSubgrid,prior,grown,current,current.capArea,v,1.,normals);
            close(subgrid(exactSubgrid)[0],grown[0],2e-15,"Spurious resolved residue from subgrid growth");
        }
        {
            SurfaceState mixed;auto& exposed=mixed.population[2];exposed.fraction=.1;exposed.count=400;exposed.supportArea=bare[2];
            MineralVector prior{},next{};prior[0]=.5;prior[2]=.4;next=prior;next[0]=.45;
            const auto surface=measure(mixed,bare,v,true,layer);
            close(inaccessibleFractions(mixed,prior,surface)[0],.4,1e-15,"Buried resolved mineral was not protected");
            grow(mixed,prior,next,surface,surface.exposed,v,1.,normals,bare,layer,ep);
            close(subgrid(mixed)[0],.05,2e-15,"Dissolution consumed buried resolved mineral");
            close(next[0]-subgrid(mixed)[0],.4,2e-15,"Buried resolved volume changed");
            auto inaccessibleAttempt=mixed;prior=next;next[0]=.3;
            bool rejected=false;try{grow(inaccessibleAttempt,prior,next,measure(mixed,bare,v,true,layer),surface.exposed,v,1.,normals,bare,layer,ep);}catch(const std::runtime_error&){rejected=true;}
            require(rejected,"Geometry accepted dissolution beyond accessible inventory");
            std::cout<<"buried_resolved_inventory_protection,1\n";
        }
        {
            SurfaceState sizes;const double r1=1e-7,r2=2e-7,advance=5e-8;
            auto volume=[&](double a){return 2*pi*a*a*a/(3*v);};
            sizes.population[2].fraction=volume(r1);sizes.population[2].count=1;
            sizes.population[3].fraction=volume(r2);sizes.population[3].count=1;
            MineralVector support=bare;support[3]=h*h;
            MineralVector prior{},next{};prior[0]=subgrid(sizes)[0];prior[2]=prior[3]=.2;
            next=prior;next[0]=volume(r1+advance)+volume(r2+advance);
            auto surface=measure(sizes,support,v);
            grow(sizes,prior,next,surface,surface.exposed,v,1.,normals,support);
            close(radius(sizes.population[2],v),r1+advance,2e-14,"Small crystal has incorrect normal growth velocity");
            close(radius(sizes.population[3],v),r2+advance,2e-14,"Large crystal has incorrect normal growth velocity");
            close(subgrid(sizes)[0],next[0],2e-15,"Mixed-radius growth created a resolved residue");
            prior=next;next[0]=volume(5e-8);surface=measure(sizes,support,v);
            grow(sizes,prior,next,surface,surface.exposed,v,1.,normals,support);
            require(sizes.population[2].fraction==0&&sizes.population[2].count==0,"Small crystal survived complete recession");
            close(radius(sizes.population[3],v),5e-8,3e-14,"Exhausted population volume was redistributed incorrectly");
            std::cout<<"population_normal_growth_and_depletion,1\n";
        }
        // Manufactured constant-affinity growth has r(t)=r0+k*Vm*t for
        // an isolated hemisphere. Integrate the production area law by RK4.
        RateParameters manufactured;manufactured.geometry.enabled=true;manufactured.geometry.separateMineralAreas=true;
        manufactured.neutral[0]=.01;manufactured.geometry.growthArea[0]=2*pi*r*r;
        manufactured.referenceMoles[0]=2*pi*r*r*r/(3*1e-4);manufactured.geometry.surfaceExponent[0]=2./3;
        MineralVector mv{},omegaSI{};mv.fill(1e-4);omegaSI[0]=std::log10(2.);
        double amount=manufactured.referenceMoles[0];const double dt=.001;
        auto rhs=[&](double n){MineralVector mineral{};mineral[0]=n;return -surfaceRatesAt(mineral,manufactured,omegaSI,-7.,-300.,mv,1.)[0];};
        for(int step=0;step<100;step++){const double k1=rhs(amount),k2=rhs(amount+.5*dt*k1),k3=rhs(amount+.5*dt*k2),k4=rhs(amount+dt*k3);amount+=dt*(k1+2*k2+2*k3+k4)/6;}
        const double finalRadius=std::cbrt(3*amount*1e-4/(2*pi));
        close(finalRadius,r+.01*1e-4*.1,1e-7,"Single-crystal radius-time law");
        std::cout<<"single_crystal_radius_relative_error,"<<std::abs(finalRadius/(r+.01*1e-4*.1)-1)<<'\n';
        Inventory cell;cell.waterEquivalent=1000;cell.aqueous.fill(1000);cell.mineral[2]=1000;
        MineralVector volumes{};volumes.fill(1e-4);SurfaceState s;
        NucleationLaw law{0,2,1e12,0,45,1e-8};law.validate();auto initial=conserved(cell);
        require(addEmbryo(cell,s,law,v,volumes,ep,0),"Embryo unexpectedly rejected");auto final=conserved(cell);
        for(int i=0;i<8;i++)close(initial[i],final[i],2e-15,"Embryo component/water balance");
        SurfaceState empty;Inventory insufficient;insufficient.aqueous.fill(0);
        require(!addEmbryo(insufficient,empty,law,v,volumes,ep,0),"Nucleation created mass without ions");
        require(law.frequency(0)==0&&law.frequency(-1)==0,"Nucleation below saturation");
        law.radius=1e-12;MineralVector si{};si[0]=1;
        {
            Inventory depleted=cell;depleted.mineral[2]=0;SurfaceState surface;
            nucleate(depleted,surface,{law},si,bare,.1,v,volumes,ep,0,2024,0);
            require(surface.births[0]==0&&surface.population[2].hazard==0,"Depleted substrate nucleated on stale area");
            require(!addEmbryo(depleted,surface,law,v,volumes,ep,0),"Embryo accepted an absent substrate");
            depleted.mineral[2]=.1/volumes[2];surface.population[14].fraction=.1;surface.population[14].count=10;
            nucleate(depleted,surface,{law},si,bare,.1,v,volumes,ep,0,2024,0);
            require(surface.births[0]==0,"Subgrid substrate nucleated on stale resolved area");
            std::cout<<"depleted_substrate_rejects_nucleation,1\n";
        }
        {
            Inventory fluid=cell;fluid.mineral[2]=0;SurfaceState surface;
            MineralVector external{};external[2]=1000;const auto substrateBefore=external;
            require(addEmbryo(fluid,surface,law,v,volumes,ep,0,bare[2],&external),"Neighbour-owned substrate could not nucleate in the liquid owner");
            require(external==substrateBefore&&fluid.mineral[2]==0,"Nucleation consumed or copied the external substrate");
            nucleate(fluid,surface,{law},si,bare,.1,v,volumes,ep,0,2024,0,{},&external);
            require(surface.births[0]>1,"External-support nucleation ignored visible area");
            MineralVector inventory{};inventory[0]=subgrid(surface)[0];const auto saved=surface;
            std::array<bool,7> support{};support[2]=true;
            grow(surface,inventory,inventory,measure(surface,bare,v),MineralVector{},v,1.,normals,bare,{},ep,&support,false);
            require(surface.population[2].count==saved.population[2].count&&surface.population[2].fraction==saved.population[2].fraction,
                "External product was promoted because its substrate is absent locally");
            auto porousExternal=surface;
            grow(porousExternal,inventory,inventory,measure(porousExternal,bare,v,true,layer),MineralVector{},v,1.,normals,bare,layer,ep,&support,true);
            require(porousExternal.population[2].retreat==0,"External product was detached because its substrate is absent locally");
            external[2]=0;const auto births=surface.births;
            nucleate(fluid,surface,{law},si,bare,.1,v,volumes,ep,0,2024,0,{},&external);
            require(surface.births==births,"Depleted external substrate continued nucleating on stale area");
            const auto detached=measure(surface,MineralVector{},v);
            require(detached.exposed[0]>0&&detached.coverage[2]==0,"Lost external support made retained impermeable products inert");
            const double independentArea=2*pi*(1-std::cos(law.contactAngle*pi/180))*std::pow(radius(surface.population[2],v),2)*surface.population[2].count;
            close(detached.capArea[0],independentArea,2e-14,"Detached product lost its independent crystal area");
            auto remaining=inventory;remaining[0]*=.5;
            grow(surface,inventory,remaining,detached,detached.exposed,v,1.,normals,{}, {},ep,nullptr,false);
            close(subgrid(surface)[0],remaining[0],2e-15,"Detached impermeable product could not dissolve");
            std::cout<<"external_substrate_guard_and_population_support,1\n";
        }
        double sum=0,sum2=0;constexpr int samples=6000;
        for(int k=0;k<samples;k++){SurfaceState t;Inventory c=cell;nucleate(c,t,{law},si,bare,.01,v,volumes,ep,0,2024,k);double n=t.births[0];sum+=n;sum2+=n*n;}
        const double mean=sum/samples,variance=sum2/samples-mean*mean;
        require(std::abs(mean-1)<.06&&std::abs(variance-1)<.1,"Poisson nucleation moments");
        std::cout<<"poisson_mean,"<<mean<<"\npoisson_variance,"<<variance<<'\n';
        SurfaceState t;Inventory c=cell;nucleate(c,t,{law},si,bare,.003,v,volumes,ep,0,99,47);
        const auto snapshot=t;const auto inventorySnapshot=c;
        nucleate(c,t,{law},si,bare,.017,v,volumes,ep,0,99,47);const auto completed=t;
        t=snapshot;c=inventorySnapshot;nucleate(c,t,{law},si,bare,.017,v,volumes,ep,0,99,47);
        require(t.births==completed.births&&t.population[2].hazard==completed.population[2].hazard,"Rejected step changed random path");
        std::cout<<"stochastic_rollback,1\n";
        // Curved-interface refinement, using independent circle quadrature
        // for volume fractions and analytic radial normals. No fitted target.
        double previous=1;
        for(int n:{16,32,64}){
            const double dh=1./n,R=.31;double perimeter=0;
            for(int iy=0;iy<n;iy++)for(int ix=0;ix<n;ix++){
                const double x=ix*dh-.5,y=iy*dh-.5;double area=0;
                for(int q=0;q<256;q++){
                    const double xx=x+(q+.5)*dh/256,z=std::sqrt(std::max(0.,R*R-xx*xx));
                    area+=std::max(0.,std::min(y+dh,z)-std::max(y,-z))*dh/256;
                }
                if(area<=dh*dh*1e-12||area>=dh*dh*(1-1e-12))continue;
                MineralVector ff{};ff[2]=area/(dh*dh)*(1-ep);auto nn=normals;nn[2]={x+.5*dh,y+.5*dh};
                perimeter+=reconstruct(ff,0,neighbours,nn,dh,dh,1,ep).internalArea[2];
            }
            const double error=std::abs(perimeter/(2*pi*R)-1);require(error<previous,"PLIC perimeter did not converge");previous=error;
            std::cout<<"circle_perimeter_error_"<<n<<','<<error<<'\n';
        }require(previous<.005,"PLIC finest-grid perimeter error");
        // Full Pitzer reaction integration with separate areas; embryo and
        // hydrate water are covered by the same conservation assertions.
        PitzerModel model(argv[1]);Speciation spec(model);Reactions reaction(model);
        {
            auto fluid=haliteSaturatedInlet(spec);const double waterScale=1000/fluid.waterEquivalent;
            for(double& x:fluid.aqueous)x*=waterScale;fluid.waterEquivalent*=waterScale;fluid.mineral[1]=.1;
            SurfaceState detached;auto& pop=detached.population[9];
            pop.fraction=fluid.mineral[1]*model.molarVolume[1];pop.count=10;pop.supportArea=h*h;pop.retreat=h;
            const auto surface=measure(detached,MineralVector{},v);
            const auto oldInventory=fluid;const auto oldBalance=conserved(fluid);
            RateParameters rp;rp.geometry.enabled=true;rp.geometry.separateMineralAreas=true;rp.neutral[1]=.02;
            rp.area[1]=rp.geometry.growthArea[1]=surface.exposed[1]/v;rp.referenceMoles[1]=.1;rp.geometry.surfaceExponent[1]=2./3;
            IntegrationOptions opts;opts.relativeTolerance=2e-9;opts.absoluteTolerance=1e-12;
            reaction.advance(fluid,.01,rp,opts);
            require(fluid.mineral[1]<oldInventory.mineral[1]&&fluid.mineral[1]>=0,"Detached impermeable product did not redissolve in Pitzer chemistry");
            const auto newBalance=conserved(fluid);for(int j=0;j<8;j++)close(oldBalance[j],newBalance[j],1e-13,"Detached impermeable product lost conserved mass");
            std::cout<<"detached_impermeable_product_redissolution,"<<fluid.mineral[1]<<'\n';
        }
        {
            auto mixed=haliteSaturatedInlet(spec);const double waterScale=1000/mixed.waterEquivalent;
            for(double& x:mixed.aqueous)x*=waterScale;mixed.waterEquivalent*=waterScale;
            mixed.aqueous[K]+=1;mixed.aqueous[Cl]+=1;
            mixed.mineral[1]=.5;const auto conservedBefore=conserved(mixed);
            RateParameters rp;rp.geometry.enabled=true;rp.geometry.separateMineralAreas=true;rp.geometry.epsilonMin=ep;
            rp.neutral[1]=.02;rp.area[1]=rp.geometry.growthArea[1]=1e4;
            rp.referenceMoles[1]=.1;rp.geometry.surfaceExponent[1]=2./3;rp.geometry.lockedMoles[1]=.4;
            auto onlyAccessible=mixed;onlyAccessible.mineral[1]-=.4;auto shifted=rp;
            shifted.geometry.lockedMoles[1]=0;shifted.geometry.inertFraction=.4*model.molarVolume[1];
            const auto evaluated=reaction.evaluate(mixed,rp),shiftedValue=reaction.evaluate(onlyAccessible,shifted);
            close(evaluated.rates[1],shiftedValue.rates[1],1e-14,"Inaccessible inventory changed active surface kinetics");
            double jac[7][7];reaction.rateJacobian(mixed,rp,evaluated,jac);
            MineralVector plus{},minus{};plus[1]=1e-6;minus[1]=-1e-6;
            const double numeric=(reaction.evaluate(Reactions::changed(mixed,plus),rp).rates[1]-reaction.evaluate(Reactions::changed(mixed,minus),rp).rates[1])/(2e-6);
            close(jac[1][1],numeric,2e-5,"Locked-inventory analytic Jacobian mismatch");
            for(int method=0;method<
#ifdef SALT_HAS_CVODE
                3
#else
                2
#endif
                ;method++){
                auto state=mixed;IntegrationOptions opts;opts.method=method;opts.relativeTolerance=2e-9;opts.absoluteTolerance=1e-14;
                reaction.advance(state,.01,rp,opts);
                require(state.mineral[1]>=.4&&state.mineral[1]<.40000001,"Integrator consumed buried inventory or failed to exhaust exposed mineral");
                const auto conservedAfter=conserved(state);for(int i=0;i<8;i++)close(conservedBefore[i],conservedAfter[i],1e-13,"Locked-inventory component/water conservation");
                require(reaction.evaluate(state,rp).rates[1]==0,"Buried inventory continued dissolving after cap depletion");
                std::cout<<"locked_inventory_integrator_"<<method<<","<<state.mineral[1]<<'\n';
            }
        }
        {
            auto residue=haliteSaturatedInlet(spec);const double waterScale=1000/residue.waterEquivalent;
            for(double& x:residue.aqueous)x*=waterScale;residue.waterEquivalent*=waterScale;
            residue.mineral[1]=.1;const auto conservedBefore=conserved(residue);
            SurfaceState particles;auto& pop=particles.population[9];pop.fraction=residue.mineral[1]*model.molarVolume[1];pop.count=10;pop.supportArea=h*h;pop.retreat=h;
            const auto surface=measure(particles,MineralVector{},v,true,layer);
            RateParameters rp;rp.geometry.enabled=true;rp.geometry.separateMineralAreas=true;rp.geometry.epsilonMin=ep;
            rp.neutral[1]=.02;rp.area[1]=rp.geometry.growthArea[1]=surface.exposed[1]/v;rp.referenceMoles[1]=.1;rp.geometry.surfaceExponent[1]=2./3;
            IntegrationOptions opts;opts.relativeTolerance=2e-9;opts.absoluteTolerance=1e-12;
            reaction.advance(residue,.01,rp,opts);
            require(residue.mineral[1]>=0&&residue.mineral[1]<.1&&residue.mineral[2]==0,"Pitzer did not redissolve retained Sylvite without Carnallite");
            const auto conservedAfter=conserved(residue);for(int i=0;i<8;i++)close(conservedBefore[i],conservedAfter[i],1e-13,"Retained-product component/water conservation");
            std::cout<<"retained_sylvite_after_redissolution,"<<residue.mineral[1]<<'\n';
        }
        auto solution=haliteSaturatedInlet(spec);const double density=spec.solve({{6,0,0,0,6,0,0}}).density;
        const double scale=.4*density/aqueousMass(solution);for(double& x:solution.aqueous)x*=scale;solution.waterEquivalent*=scale;solution.mineral[2]=.6/model.molarVolume[2];
        RateParameters rates;rates.geometry.enabled=true;rates.geometry.separateMineralAreas=true;rates.geometry.bulkVolume=1;rates.geometry.epsilonMin=ep;
        rates.neutral[2]=.02;rates.area[2]=1/h;rates.referenceMoles[2]=solution.mineral[2];rates.geometry.growthArea[2]=1/h;
        initial=conserved(solution);IntegrationOptions options;options.relativeTolerance=2e-9;options.absoluteTolerance=1e-12;
        reaction.advance(solution,.01,rates,options);final=conserved(solution);double error=0;
        require(solution.mineral[0]==0&&solution.mineral[1]==0,"Implicit solver seeded an inactive mineral");
        for(int i=0;i<8;i++)error=std::max(error,std::abs(initial[i]-final[i])/std::max(std::abs(initial[i]),1.));require(error<1e-13,"Pitzer interface mass balance");
        std::cout<<"pitzer_interface_balance,"<<error<<'\n';
        solution.mineral[0]=.1;rates.neutral[0]=.02;rates.area[0]=rates.geometry.growthArea[0]=1e4;rates.referenceMoles[0]=.1;rates.geometry.surfaceExponent[0]=2./3;
        const auto evaluated=reaction.evaluate(solution,rates);double jac[7][7];reaction.rateJacobian(solution,rates,evaluated,jac);double derivativeError=0;
        for(int j:{0,2}){
            const double step=1e-6*std::max(solution.mineral[j],1.);MineralVector plus{},minus{};plus[j]=step;minus[j]=-step;
            const auto a=reaction.evaluate(Reactions::changed(solution,plus),rates).rates,b=reaction.evaluate(Reactions::changed(solution,minus),rates).rates;
            for(int k:{0,2}){const double numeric=(a[k]-b[k])/(2*step);derivativeError=std::max(derivativeError,std::abs(jac[k][j]-numeric)/std::max({std::abs(jac[k][j]),std::abs(numeric),1e-10}));}
        }
        require(derivativeError<3e-4,"Separate-area analytic Jacobian mismatch");
        std::cout<<"separate_area_jacobian_error,"<<derivativeError<<"\nall_passed,1\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
