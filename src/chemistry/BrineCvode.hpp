// SPDX-License-Identifier: GPL-3.0-or-later
// Included after Reactions is defined. Only the ODE time integrator is external;
// Pitzer speciation, TST rates, geometry and conservation remain native code.
#pragma once
namespace saltchem {
inline IntegrationStats Reactions::advanceCvode(Inventory& state,double interval,const RateParameters& p,const IntegrationOptions& options)const {
    IntegrationStats stats;if(interval==0)return stats;const Inventory intervalStart=state;
    struct Workspace {
        SUNContext context=nullptr;N_Vector y=nullptr,constraints=nullptr;SUNMatrix matrix=nullptr;
        SUNLinearSolver solver=nullptr;void* memory=nullptr;
        const Reactions* reaction=nullptr;const RateParameters* parameters=nullptr;
        const IntegrationOptions* options=nullptr;IntegrationStats* stats=nullptr;
        Inventory base;RateParameters segmentParameters;AqueousState seed;std::string diagnostic;
        static void checked(int code){if(code<0)throw std::runtime_error("CVODE setup error "+std::to_string(code));}
        static MineralVector vector(N_Vector v){MineralVector a{};for(int j=0;j<nMinerals;j++)a[j]=NV_Ith_S(v,j);return a;}
        ReactionOutput eval(const MineralVector& proposed){
            auto extent=proposed;auto s=Reactions::changed(base,extent);
            // Keep the finite-precision zero-phase boundary inside the local
            // error allowance. The paired dissolved transfer is corrected by
            // changing extent, never by independently clipping inventories.
            for(int j=0;j<nMinerals;j++)if(s.mineral[j]<0&&s.mineral[j]>=-options->absoluteTolerance)extent[j]=base.mineral[j];
            s=Reactions::changed(base,extent);
            if(!Reactions::valid(s))throw std::domain_error("CVODE trial violates positivity");
            auto out=reaction->evaluate(s,*parameters,&seed,false,reaction->incrementedPorosity(base,extent,*parameters));
            seed=out.aqueous;stats->rateEvaluations++;stats->speciationIterations+=out.aqueous.iterations;return out;
        }
        static int rhs(sunrealtype,N_Vector y,N_Vector dy,void* data){
            auto& w=*static_cast<Workspace*>(data);
            try{auto out=w.eval(vector(y));for(int j=0;j<nMinerals;j++)NV_Ith_S(dy,j)=out.rates[j];NV_Ith_S(dy,nMinerals)=0;return 0;}
            catch(const std::exception& e){w.diagnostic=e.what();return 1;}
            catch(...){return -1;}
        }
        static int weights(N_Vector y,N_Vector weights,void* data){
            auto& w=*static_cast<Workspace*>(data);const auto s=Reactions::changed(w.base,vector(y));
            const auto& o=*w.options;
            // Conservative extent weights bound errors in every aqueous total,
            // solid inventory and water, including nearly closed micropores.
            // The sqrt(dimension) factor converts the RMS to a maximum norm.
            for(int j=0;j<nMinerals;j++){
                double budget=o.absoluteTolerance+o.relativeTolerance*std::abs(s.mineral[j]);
                for(int i=0;i<nComponents;i++)if(stoichiometry[j][i]){
                    double rowSum=0;for(int k=0;k<nMinerals;k++)rowSum+=std::abs(stoichiometry[k][i]);
                    budget=std::min(budget,(o.absoluteTolerance+o.relativeTolerance*std::abs(s.aqueous[i]))/rowSum);
                }
                if(hydration[j]){double sum=0;for(double h:hydration)sum+=h*waterMolarMass;budget=std::min(budget,(waterMolarMass*o.absoluteTolerance+o.relativeTolerance*std::abs(s.waterEquivalent))/sum);}
                if(!(budget>0)||!std::isfinite(budget))return -1;
                NV_Ith_S(weights,j)=std::sqrt(double(nMinerals+1))/budget;
            }NV_Ith_S(weights,nMinerals)=1;return 0;
        }
        static int jacobian(sunrealtype,N_Vector y,N_Vector fy,SUNMatrix jac,void* data,N_Vector,N_Vector,N_Vector){
            auto& w=*static_cast<Workspace*>(data);const auto x=vector(y),initial=vector(fy);const auto s=Reactions::changed(w.base,x);SUNMatZero(jac);
            if(w.options->analyticJacobian){
                try{
                    const auto value=w.eval(x);double derivative[nMinerals][nMinerals];
                    w.reaction->rateJacobian(s,*w.parameters,value,derivative,w.reaction->incrementedPorosity(w.base,x,*w.parameters));
                    for(int i=0;i<nMinerals;i++)for(int j=0;j<nMinerals;j++)SM_ELEMENT_D(jac,i,j)=derivative[i][j];
                    w.stats->analyticJacobians++;return 0;
                }catch(const std::exception& e){w.diagnostic=e.what();w.stats->jacobianFallbacks++;}
            }
            for(int j=0;j<nMinerals;j++){
                for(int i=0;i<nMinerals;i++)SM_ELEMENT_D(jac,i,j)=0;
                const auto& p=*w.parameters;
                if((!p.geometry.enabled&&p.area[j]==0)||(p.geometry.enabled&&p.geometry.weights[j]==0))continue;
                double delta=1e-6*std::max(std::min(s.mineral[j],s.waterEquivalent),1e-6*s.waterEquivalent);bool done=false;
                for(int attempt=0;attempt<20;attempt++,delta*=.1){
                    auto plus=x,minus=x;plus[j]+=delta;minus[j]-=delta;
                    const bool vp=w.reaction->admissible(Reactions::changed(w.base,plus),p),vm=w.reaction->admissible(Reactions::changed(w.base,minus),p);
                    if(!(vp||vm))continue;
                    try{const auto a=vp?w.eval(plus).rates:initial,b=vm?w.eval(minus).rates:initial;
                        for(int i=0;i<nMinerals;i++)SM_ELEMENT_D(jac,i,j)=(a[i]-b[i])/((vp&&vm?2:1)*delta);done=true;break;
                    }catch(const std::exception& e){w.diagnostic=e.what();}
                }
                if(!done&&s.mineral[j]>0)return 1;
            }return 0;
        }
        static int roots(sunrealtype,N_Vector y,sunrealtype* out,void* data){
            auto& w=*static_cast<Workspace*>(data);
            for(int j=0;j<nMinerals;j++)out[j]=w.base.mineral[j]-NV_Ith_S(y,j)-w.options->absoluteTolerance;
            return 0;
        }
        Workspace(){
            checked(SUNContext_Create(SUN_COMM_NULL,&context));
            checked(SUNContext_ClearErrHandlers(context));
            y=N_VNew_Serial(nMinerals+1,context);constraints=N_VNew_Serial(nMinerals+1,context);matrix=SUNDenseMatrix(nMinerals+1,nMinerals+1,context);
            solver=SUNLinSol_Dense(y,matrix,context);memory=CVodeCreate(CV_BDF,context);
            if(!y||!constraints||!matrix||!solver||!memory)throw std::bad_alloc();
            N_VConst(0,y);checked(CVodeInit(memory,rhs,0,y));checked(CVodeSetUserData(memory,this));
            checked(CVodeWFtolerances(memory,weights));checked(CVodeSetLinearSolver(memory,solver,matrix));checked(CVodeSetJacFn(memory,jacobian));
            checked(CVodeRootInit(memory,nMinerals,roots));int direction[nMinerals];std::fill_n(direction,nMinerals,-1);checked(CVodeSetRootDirection(memory,direction));
            checked(CVodeSetMaxErrTestFails(memory,20));checked(CVodeSetMaxConvFails(memory,20));checked(CVodeSetMaxNonlinIters(memory,6));
        }
        ~Workspace(){if(memory)CVodeFree(&memory);if(solver)SUNLinSolFree(solver);if(matrix)SUNMatDestroy(matrix);if(constraints)N_VDestroy(constraints);if(y)N_VDestroy(y);if(context)SUNContext_Free(&context);}
        void counters(){long int n=0;CVodeGetNumErrTestFails(memory,&n);stats->rejected+=int(n);CVodeGetNumNonlinSolvConvFails(memory,&n);stats->rejected+=int(n);CVodeGetNumNonlinSolvIters(memory,&n);stats->newtonIterations+=int(n);}
    };
    // Allocation and matrix storage are reused by each OpenMP worker. ReInit
    // discards time history when the external transport state changes.
    static thread_local Workspace w;w.reaction=this;w.segmentParameters=p;w.parameters=&w.segmentParameters;w.options=&options;w.stats=&stats;w.base=state;w.diagnostic.clear();
    double porosity=incrementedPorosity(state,{},p);
    auto initial=evaluate(state,p,nullptr,false,porosity);w.seed=initial.aqueous;stats.rateEvaluations++;stats.speciationIterations+=initial.aqueous.iterations;
    if(std::all_of(initial.rates.begin(),initial.rates.end(),[](double r){return r==0;}))return stats;
    double time=0;int calls=0;MineralVector previous{};
    auto restart=[&](){
        w.base=state;w.segmentParameters.geometry.initialPorosity=porosity;previous.fill(0);N_VConst(0,w.y);
        // Separate surface operators can permit only one reaction direction.
        // Enforce that sign on the integrated extent as well as on the RHS:
        // a dense BDF solve can otherwise create a tiny opposite-direction
        // extent in a saturated phase coupled to another dissolving mineral.
        // The passive zero variable keeps the constraint vector nonempty.
        for(int j=0;j<nMinerals;j++){
            const bool separate=p.geometry.separateMineralAreas;
            NV_Ith_S(w.constraints,j)=separate&&p.geometry.growthArea[j]==0?1:
                (state.mineral[j]==0||(separate&&p.area[j]==0)?-1:0);
        }
        NV_Ith_S(w.constraints,nMinerals)=1;
        Workspace::checked(CVodeReInit(w.memory,time,w.y));Workspace::checked(CVodeSetConstraints(w.memory,w.constraints));
        Workspace::checked(CVodeSetStopTime(w.memory,interval));Workspace::checked(CVodeSetMaxNumSteps(w.memory,options.maximumAttempts));
    };
    restart();
    auto finishWithSingleStep=[&](){
        // Recompute the whole local interval if the multistep path stalls.
        // A partially advanced stiff state may carry the same history-dependent
        // difficulty into another method. The external conserved state is only
        // committed after one method has completed the entire interval.
        w.counters();state=intervalStart;const auto tail=advanceNormalized(state,interval,p,options);
        stats.integratorFallbacks++;stats.accepted+=tail.accepted;stats.rejected+=tail.rejected;
        stats.rateEvaluations+=tail.rateEvaluations;stats.speciationIterations+=tail.speciationIterations;
        stats.newtonIterations+=tail.newtonIterations;stats.depletionEvents=tail.depletionEvents;
        stats.dissolved=tail.dissolved;stats.precipitated=tail.precipitated;
        return stats;
    };
    while(time<interval*(1-1e-13)){
        // Preserve a single-step method for rare intervals where the BDF history
        // becomes inefficient during interface evolution. The alternative uses
        // the same external interval and error tolerances.
        if(calls>=options.maximumAttempts)return finishWithSingleStep();
        ++calls;
        const double lastTime=time;
        const int code=CVode(w.memory,interval,w.y,&time,CV_ONE_STEP);
        if(code<0){
            // A changing interface can invalidate the multistep history. Rebuild
            // it at the last accepted, conserved state, with the same tolerances.
            if((code==CV_CONV_FAILURE||code==CV_ERR_FAILURE)&&stats.historyRestarts<6){stats.historyRestarts++;w.counters();time=lastTime;w.seed=evaluate(state,p,nullptr,false,porosity).aqueous;restart();continue;}
            if(code==CV_CONV_FAILURE||code==CV_ERR_FAILURE){time=lastTime;return finishWithSingleStep();}
            std::ostringstream error;error<<std::setprecision(15)<<"CVODE code="<<code<<" t="<<time<<" steps="<<stats.accepted<<" restarts="<<stats.historyRestarts<<" detail="<<w.diagnostic;throw std::runtime_error(error.str());
        }
        auto extent=Workspace::vector(w.y);state=changed(w.base,extent);bool repaired=false;double correction=0;
        // Newton's converged iterate may cross an exactly absent solid by less
        // than the absolute error allowance. Correct its reaction extent, so
        // the matching aqueous/water transfer is corrected conservatively too.
        for(int j=0;j<nMinerals;j++)if(state.mineral[j]<0&&state.mineral[j]>=-options.absoluteTolerance){correction=std::max(correction,-state.mineral[j]);extent[j]=w.base.mineral[j];repaired=true;}
        if(repaired)state=changed(w.base,extent);
        porosity=incrementedPorosity(w.base,extent,*w.parameters);
        if(!admissible(state,p)){std::ostringstream error;error<<std::setprecision(17)<<"CVODE accepted an inadmissible state t="<<time<<" epsilon="<<poreFraction(state,p.geometry,model.molarVolume)<<" water="<<state.waterEquivalent<<" aqueous=";for(double v:state.aqueous)error<<v<<',';error<<" minerals=";for(double v:state.mineral)error<<v<<',';throw std::runtime_error(error.str());}
        for(int j=0;j<nMinerals;j++){const double change=extent[j]-previous[j];stats.dissolved[j]+=std::max(change,0.);stats.precipitated[j]+=std::max(-change,0.);}
        previous=extent;stats.accepted++;
        if(code==CV_ROOT_RETURN){
            int roots[nMinerals];Workspace::checked(CVodeGetRootInfo(w.memory,roots));MineralVector remainder{};
            for(int j=0;j<nMinerals;j++)if(roots[j]<0){remainder[j]=state.mineral[j];stats.dissolved[j]+=remainder[j];stats.depletionEvents++;}
            auto eventParameters=p;eventParameters.geometry.initialPorosity=porosity;
            porosity=incrementedPorosity(state,remainder,eventParameters);
            state=changed(state,remainder);w.counters();restart();
        }else if(correction>.01*options.absoluteTolerance){w.counters();restart();}
    }
    w.counters();return stats;
}
}
