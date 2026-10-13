// SPDX-License-Identifier: GPL-3.0-or-later
#include "BrineReactions.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
using namespace saltchem;
int main(int argc,char** argv){try{
    const std::string db=argc>1?argv[1]:"src/chemistry/pitzer.dat";
    PitzerModel model(db);Speciation spec(model);Reactions reaction(model);
    const Totals totals{{4,.2,.1,.5,5.328,.035,.002}};MineralVector minerals;minerals.fill(.2);
    RateParameters p;p.neutral={{.02,.002,.02,1e-5,1e-4,.02,1e-5}};p.area.fill(.1);p.referenceMoles=minerals;
    double worst=0;int count=0;
    for(auto tp:std::initializer_list<std::array<double,2>>{{25,1},{60,100},{150,500}}){
        model.setConditions(tp[0]+273.15,tp[1]);
        for(bool geometry:{false,true}){
            auto v=spec.inventory(totals,1,minerals);p.geometry.enabled=geometry;
            p.geometry.bulkVolume=.002;p.geometry.inertFraction=.6;p.geometry.neighbours={{.2,.7,.25,.65}};
            const auto base=reaction.evaluate(v,p);double jac[7][7];reaction.rateJacobian(v,p,base,jac);
            for(int j=0;j<7;j++){
                MineralVector shift{};const double delta=3e-6;shift[j]=delta;
                const auto a=reaction.evaluate(Reactions::changed(v,shift),p,&base.aqueous).rates;shift[j]=-delta;
                const auto b=reaction.evaluate(Reactions::changed(v,shift),p,&base.aqueous).rates;
                for(int i=0;i<7;i++){
                    double rowScale=0;for(int k=0;k<7;k++)rowScale=std::max(rowScale,std::abs(jac[i][k]));
                    const double fd=(a[i]-b[i])/(2*delta),error=std::abs(jac[i][j]-fd)/(1e-8+rowScale);
                    worst=std::max(worst,error);if(error>3e-5)std::cerr<<"derivative "<<tp[0]<<' '<<geometry<<' '<<i<<' '<<j<<' '<<jac[i][j]<<' '<<fd<<'\n';count++;
                }
            }
        }
    }
    p.geometry.enabled=false;model.setConditions(298.15,1);auto initial=spec.inventory(totals,1,minerals);
    double times[2]={},stateError=0;Inventory answers[2];long long evaluations[2]={},jacobians=0,fallbacks=0;
    for(int mode=0;mode<2;mode++){
        IntegrationOptions options;options.analyticJacobian=mode==1;const auto start=std::chrono::steady_clock::now();
        for(int n=0;n<60;n++){auto v=initial;const auto st=reaction.advance(v,.5,p,options);answers[mode]=v;evaluations[mode]+=st.rateEvaluations;if(mode){jacobians+=st.analyticJacobians;fallbacks+=st.jacobianFallbacks;}}
        times[mode]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    for(int j=0;j<7;j++){stateError=std::max(stateError,std::abs(answers[0].aqueous[j]-answers[1].aqueous[j]));stateError=std::max(stateError,std::abs(answers[0].mineral[j]-answers[1].mineral[j]));}
    std::cout<<std::setprecision(14)<<"derivative_entries="<<count<<" max_scaled_error="<<worst<<" state_difference="<<stateError<<" finite_difference_s="<<times[0]<<" analytic_s="<<times[1]<<" FD_RHS="<<evaluations[0]<<" AD_RHS="<<evaluations[1]<<" analytic_jacobians="<<jacobians<<" fallbacks="<<fallbacks<<'\n';
    return worst>3e-5||stateError>1e-9||jacobians==0||fallbacks!=0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
