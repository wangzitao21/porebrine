// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
namespace saltchem {
// Implicit differentiation of the converged activity/charge/solvent equations.
// No differentiation through a finite number of speciation iterations.
inline void Reactions::rateJacobian(const Inventory& state,const RateParameters& p,const ReactionOutput& value,double result[nMinerals][nMinerals],double porosity)const {
    constexpr int nx=nSpecies+4, nd=nx+nMinerals;
    using D=Dual<nd>;using R=Dual<nMinerals>;
    const auto& s=value.aqueous;
    std::array<bool,nMinerals> feasibleDirection;
    for(int j=0;j<nMinerals;j++){
        feasibleDirection[j]=true;
        if(p.geometry.separateMineralAreas&&p.area[j]==0&&p.geometry.growthArea[j]==0)feasibleDirection[j]=false;
        if(state.mineral[j]<=p.geometry.lockedMoles[j])for(int i=0;i<nComponents;i++)if(stoichiometry[j][i]&&state.aqueous[i]==0)feasibleDirection[j]=false;
    }
    std::array<D,nx> x;
    for(int i=0;i<nSpecies;i++)x[i]=D::variable(s.logGamma[i],i);
    x[nSpecies]=D::variable(std::log(s.waterActivity),nSpecies);
    x[nSpecies+1]=D::variable(-s.pH,nSpecies+1);
    x[nSpecies+2]=D::variable(s.ionicStrength,nSpecies+2);
    x[nSpecies+3]=D::variable(std::log(value.freeWaterKg),nSpecies+3);
    std::array<D,nComponents> amount{},totals{};D water=state.waterEquivalent;
    for(int i=0;i<nComponents;i++)amount[i]=state.aqueous[i];
    for(int j=0;j<nMinerals;j++){
        const D extent=feasibleDirection[j]?D::variable(0,nx+j):D(0);
        for(int i=0;i<nComponents;i++)amount[i]+=stoichiometry[j][i]*extent;
        water+=hydration[j]*waterMolarMass*extent;
    }
    const D freeWater=exp(x[nSpecies+3]);
    for(int i=0;i<nComponents;i++)totals[i]=amount[i]/freeWater;
    std::array<D,nSpecies> gamma;for(int i=0;i<nSpecies;i++)gamma[i]=x[i];
    const auto m=Speciation::distribute(model,totals,gamma,x[nSpecies],x[nSpecies+2],x[nSpecies+1]);
    // Pitzer depends on 14 species and aw, not on all 25 residual/input
    // directions independently. Differentiate that smaller graph once and
    // apply the chain rule to the distribution derivatives.
    std::array<int,nSpecies> inputs{};int count=0;
    for(int i=0;i<nSpecies;i++)if(std::any_of(m[i].derivative.begin(),m[i].derivative.end(),[](double d){return d!=0;}))inputs[count++]=i;
    const D solventAw=exp(x[nSpecies]);std::array<D,nSpecies> computed;D lnAw,I;
    auto activityDerivative=[&](auto tag){
        using A=decltype(tag);std::array<A,nSpecies> localM,localGamma;
        for(int i=0;i<nSpecies;i++)localM[i]=m[i].value;
        for(int i=0;i<count;i++)localM[inputs[i]]=A::variable(m[inputs[i]].value,i);
        A localAw,localI;
        model.activities(localM,localGamma,localAw,localI,A::variable(solventAw.value,count));
        auto compose=[&](const A& q){D out(q.value);for(int k=0;k<nd;k++){
            double d=q.derivative[count]*solventAw.derivative[k];for(int i=0;i<count;i++)d+=q.derivative[i]*m[inputs[i]].derivative[k];out.derivative[k]=d;
        }return out;};
        for(int i=0;i<nSpecies;i++)computed[i]=compose(localGamma[i]);lnAw=compose(localAw);I=compose(localI);
    };
    // Exactly absent components with no admissible reaction direction have
    // zero derivatives. The smaller graph uses no concentration threshold.
    if(count<8)activityDerivative(Dual<8>{});else activityDerivative(Dual<nSpecies+1>{});
    std::array<D,nx> residual;
    for(int i=0;i<nSpecies;i++)residual[i]=x[i]-computed[i];
    residual[nSpecies]=x[nSpecies]-lnAw;
    residual[nSpecies+1]=totals[Na]+totals[K]+2*totals[Ca]+2*totals[Mg]-totals[Cl]-2*totals[S]-2*totals[C]+m[H]-m[OH]+m[Bicarbonate]+2*m[CO2]+m[Bisulfate]-m[MgOH];
    residual[nSpecies+2]=(x[nSpecies+2]-I)/(1+s.ionicStrength);
    residual[nSpecies+3]=(freeWater*(1+waterMolarMass*(m[OH]+m[MgOH]-m[CO2]))-water)/state.waterEquivalent;
    double a[nx][nx+nMinerals];
    for(int i=0;i<nx;i++){
        double scale=0;for(int k=0;k<nx;k++)scale=std::max(scale,std::abs(residual[i].derivative[k]));
        if(!(scale>0)||!std::isfinite(scale))throw std::runtime_error("Singular speciation derivative row");
        for(int k=0;k<nx;k++)a[i][k]=residual[i].derivative[k]/scale;
        for(int j=0;j<nMinerals;j++)a[i][nx+j]=-residual[i].derivative[nx+j]/scale;
    }
    for(int k=0;k<nx;k++){
        int pivot=k;for(int i=k+1;i<nx;i++)if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
        if(std::abs(a[pivot][k])<1e-18||!std::isfinite(a[pivot][k]))throw std::runtime_error("Singular speciation derivative matrix");
        if(pivot!=k)for(int j=k;j<nx+nMinerals;j++)std::swap(a[k][j],a[pivot][j]);
        for(int i=k+1;i<nx;i++){
            const double factor=a[i][k]/a[k][k];
            for(int j=k+1;j<nx+nMinerals;j++)a[i][j]-=factor*a[k][j];
        }
    }
    std::array<R,nx> totalX;for(int i=0;i<nx;i++)totalX[i]=x[i].value;
    for(int j=0;j<nMinerals;j++)for(int i=nx-1;i>=0;i--){
        double v=a[i][nx+j];for(int k=i+1;k<nx;k++)v-=a[i][k]*totalX[k].derivative[j];
        totalX[i].derivative[j]=v/a[i][i];
    }
    std::array<R,nComponents> totalAmount{},totalMolality{};
    std::array<R,nMinerals> mineral{},si{},omega{};
    R ep=p.geometry.enabled?(std::isfinite(porosity)?porosity:poreFraction(state,p.geometry,model.molarVolume)):1;
    for(int i=0;i<nComponents;i++){
        totalAmount[i]=state.aqueous[i];
        for(int j=0;j<nMinerals;j++)if(feasibleDirection[j])totalAmount[i].derivative[j]=stoichiometry[j][i];
        totalMolality[i]=totalAmount[i]/exp(totalX[nSpecies+3]);
    }
    for(int j=0;j<nMinerals;j++){
        mineral[j]=state.mineral[j];
        if(feasibleDirection[j]){mineral[j].derivative[j]=-1;if(p.geometry.enabled)ep.derivative[j]=model.molarVolume[j]/p.geometry.bulkVolume;}
    }
    std::array<R,nSpecies> totalGamma;for(int i=0;i<nSpecies;i++)totalGamma[i]=totalX[i];
    const auto totalSpecies=Speciation::distribute(model,totalMolality,totalGamma,totalX[nSpecies],totalX[nSpecies+2],totalX[nSpecies+1]);
    const auto k=model.phaseLogK(totalX[nSpecies+2],exp(totalX[nSpecies]));
    for(int j=0;j<nMinerals;j++){
        si[j]=s.saturationIndex[j];
        omega[j]=pow10(-k[j]+hydration[j]*totalX[nSpecies]/ln10);
        for(int i=0;i<nComponents;i++)if(stoichiometry[j][i])omega[j]*=pow(totalSpecies[i]*pow10(totalGamma[i]),stoichiometry[j][i]);
    }
    const R co2Activity=totalSpecies[CO2]*pow10(totalGamma[CO2]);
    const auto rates=surfaceRatesAt(mineral,p,si,totalX[nSpecies+1],R(-300),model.molarVolume,ep,static_cast<R*>(nullptr),static_cast<R*>(nullptr),&omega,&co2Activity);
    for(int i=0;i<nMinerals;i++)for(int j=0;j<nMinerals;j++){
        result[i][j]=rates[i].derivative[j];
        if(!std::isfinite(result[i][j]))throw std::runtime_error("Nonfinite kinetic Jacobian");
    }
}
}
