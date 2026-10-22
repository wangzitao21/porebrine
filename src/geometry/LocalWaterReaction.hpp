// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "DualSurfaceChemistry.hpp"

namespace saltdual {
struct WaterTask {std::vector<int> regions,links;};
inline std::vector<WaterTask> waterTasks(const Geometry& g,const ReactionMap& map){
    std::vector<int> parent(g.waterRegions.size());std::iota(parent.begin(),parent.end(),0);
    auto root=[&](int a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;};
    std::vector<unsigned char> needed(parent.size());for(int r:map.active)needed[r]=1;
    for(const auto& edge:g.waterLinks)if(resolvedWater(g,edge.a)&&resolvedWater(g,edge.b)){parent[root(edge.a)]=root(edge.b);needed[edge.a]=needed[edge.b]=1;}
    std::map<int,int> ids;std::vector<WaterTask> tasks;
    for(int r=0;r<int(parent.size());r++)if(needed[r]){auto inserted=ids.emplace(root(r),tasks.size());if(inserted.second)tasks.emplace_back();tasks[inserted.first->second].regions.push_back(r);}
    for(int k=0;k<int(g.waterLinks.size());k++)if(resolvedWater(g,g.waterLinks[k].a)&&resolvedWater(g,g.waterLinks[k].b))tasks[ids.at(root(g.waterLinks[k].a))].links.push_back(k);
    return tasks;
}
inline bool solveDense(std::vector<double> a,std::vector<double> rhs,std::vector<double>& x){
    int n=rhs.size();x.assign(n,0);
    for(int k=0;k<n;k++){
        int pivot=k;for(int i=k+1;i<n;i++)if(std::abs(a[i*n+k])>std::abs(a[pivot*n+k]))pivot=i;
        if(!std::isfinite(a[pivot*n+k])||std::abs(a[pivot*n+k])<1e-24)return false;
        if(pivot!=k){for(int j=k;j<n;j++)std::swap(a[k*n+j],a[pivot*n+j]);std::swap(rhs[k],rhs[pivot]);}
        for(int i=k+1;i<n;i++){double factor=a[i*n+k]/a[k*n+k];for(int j=k+1;j<n;j++)a[i*n+j]-=factor*a[k*n+j];rhs[i]-=factor*rhs[k];}
    }
    for(int i=n-1;i>=0;i--){double aij=rhs[i];for(int j=i+1;j<n;j++)aij-=a[i*n+j]*x[j];x[i]=aij/a[i*n+i];}return true;
}

// Reaction changes both dissolved mass and pore volume. Connected retained
// water therefore draws in or expels solution, in addition to diffusion.
// Eliminate its quasi-steady pressure locally with the bulk nodes as pressure
// references. Upwind transport is triangular in pressure order; every amount
// leaving one reservoir enters its neighbour, including hydrate water.
inline void redistributeReactionWater(const Geometry& exchange,const std::vector<unsigned char>& retained,
    const std::vector<double>& targetMass,std::vector<saltchem::Inventory>& state){
    const int n=state.size();std::vector<int> index(n,-1);int count=0;
    auto constrained=retained;auto target=targetMass;
    // A sealed connected pore still exchanges water internally. Fix only a
    // pressure gauge; distribute its finite mass according to the common
    // density-volume constraint, without adding an external water source.
    if(std::all_of(retained.begin(),retained.end(),[](unsigned char v){return v!=0;})){
        long double available=0,desired=0;int gauge=0;
        for(int i=0;i<n;i++){available+=saltchem::aqueousMass(state[i]);desired+=target[i];if(target[i]>target[gauge])gauge=i;}
        if(!(available>0&&desired>0))throw StepFailure("Closed local water has exhausted its finite inventory");
        for(double& value:target)value*=double(available/desired);
        constrained[gauge]=0;
    }
    for(int i=0;i<n;i++)if(constrained[i])index[i]=count++;
    if(count==0)return;
    std::vector<double> matrix(count*count),rhs(count),pressure(n),solution;
    for(int i=0;i<n;i++)if(constrained[i])rhs[index[i]]=saltchem::aqueousMass(state[i])-target[i];
    auto conductance=[](const WaterLink& edge){return edge.hydraulicConductance>0?edge.hydraulicConductance:edge.conductance;};
    double scale=0;for(const auto& edge:exchange.waterLinks)scale=std::max(scale,conductance(edge));
    if(!(scale>0))return;
    for(const auto& e:exchange.waterLinks){double k=conductance(e)/scale;int a=index[e.a],b=index[e.b];
        if(a>=0){matrix[a*count+a]+=k;if(b>=0)matrix[a*count+b]-=k;}
        if(b>=0){matrix[b*count+b]+=k;if(a>=0)matrix[b*count+a]-=k;}
    }
    if(!solveDense(matrix,rhs,solution))throw StepFailure("Singular local displacement-water network");
    for(int i=0;i<n;i++)if(constrained[i])pressure[i]=solution[index[i]];
    struct Flux {int to;double mass;};std::vector<std::vector<Flux>> outgoing(n);
    std::vector<double> mass(n),out(n);std::vector<WaterStock> stocks(n);
    for(int i=0;i<n;i++){stocks[i]=waterStock(state[i]);mass[i]=waterMass(stocks[i]);}
    for(const auto& e:exchange.waterLinks){double q=conductance(e)/scale*(pressure[e.a]-pressure[e.b]);int a=q>=0?e.a:e.b,b=q>=0?e.b:e.a;q=std::abs(q);
        if(q>0){outgoing[a].push_back({b,q});out[a]+=q;mass[a]-=q;mass[b]+=q;}}
    std::vector<int> order(n);std::iota(order.begin(),order.end(),0);std::sort(order.begin(),order.end(),[&](int a,int b){return pressure[a]>pressure[b];});
    for(int i:order){
        if(!(mass[i]>0&&mass[i]+out[i]>0))throw StepFailure("Local displacement flow exhausted its supplying water");
        auto fraction=stocks[i];for(double& a:fraction)a/=mass[i]+out[i];
        for(const auto& flux:outgoing[i])addWater(stocks[flux.to],fraction,flux.mass);
        for(int j=0;j<8;j++)stocks[i][j]=fraction[j]*mass[i];setWater(state[i],stocks[i]);
    }
}

// Diffusion is eliminated with its common conservative M-matrix. Newton
// solves only the active mineral extents. This couples water exchange and
// chemistry implicitly, including the vanishing-water/thin-film limit.
// Ordinary reservoirs keep the existing fast CVODE path.
template<class Worker>
std::vector<RegionReaction> reactWaterTask(const WaterTask& task,const Geometry& g,const ReactionMap& map,
    const std::vector<saltchem::Inventory>& cells,const saltchem::RateParameters& base,const saltchem::MineralVector& volumes,
    double diffusivity,double dt,const saltchem::IntegrationOptions& options,Worker worker){
    std::vector<RegionReaction> result;
    const bool hasRetained=std::any_of(task.regions.begin(),task.regions.end(),[&](int r){return !g.waterRegions[r].bulk;});
    if(task.links.empty()||(diffusivity==0&&!hasRetained)){for(int r:task.regions)result.push_back(reactRegion(r,g,map,cells,worker(g.waterRegions[r].parent).reactions,base,volumes,dt,options));return result;}
    const int n=task.regions.size();std::vector<RegionProblem> problems;std::vector<int> local(g.waterRegions.size(),-1);
    Geometry exchange;exchange.nx=n;exchange.ny=exchange.ratio=1;exchange.spacing=g.spacing;exchange.depth=g.depth;exchange.epsilon=0;
    std::vector<double> initialVolume,initialDensity;std::vector<unsigned char> retained(n);std::vector<saltchem::IntegrationStats> stats(n);
    std::vector<std::pair<int,int>> variables;
    for(int i=0;i<n;i++){
        int r=task.regions[i];local[r]=i;problems.push_back(regionProblem(r,g,map,cells,base,volumes));
        double v=g.waterRegions[r].fraction+(g.waterRegions[r].bulk?g.epsilon:0);initialVolume.push_back(v);
        initialDensity.push_back(saltchem::aqueousMass(problems.back().initial)/v);retained[i]=!g.waterRegions[r].bulk;
        exchange.waterRegions.push_back({i,true,true,v,{},{}});
        for(int m=0;m<3;m++)if(base.neutral[m]>0&&(problems.back().rates.area[m]>0||problems.back().rates.geometry.growthArea[m]>0))variables.push_back({i,m});
    }
    for(int k:task.links){const auto& link=g.waterLinks[k];exchange.waterLinks.push_back({local[link.a],local[link.b],link.conductance,link.hydraulicConductance});}
    int dimension=variables.size();std::vector<double> x(dimension),f(dimension),jac(dimension*dimension),direction;
    // Separate a retained pool's net volume change from exchange between
    // growing minerals. This is an invertible Newton coordinate transform;
    // the nonlinear equations, rates and convergence tolerances are unchanged.
    std::vector<int> pivot(n,-1),anchor(dimension);
    for(int k=0;k<dimension;k++){auto [i,m]=variables[k];if(!g.waterRegions[task.regions[i]].bulk&&(pivot[i]<0||volumes[m]>volumes[variables[pivot[i]].second]))pivot[i]=k;}
    for(int k=0;k<dimension;k++)anchor[k]=pivot[variables[k].first]<0?k:pivot[variables[k].first];
    auto toPhysical=[&](const std::vector<double>& y){std::vector<long double> wide(y.begin(),y.end());for(int k=0;k<dimension;k++)if(anchor[k]!=k)wide[anchor[k]]-=static_cast<long double>(volumes[variables[k].second])/volumes[variables[anchor[k]].second]*y[k];return std::vector<double>(wide.begin(),wide.end());};
    auto toBasis=[&](const std::vector<double>& values){std::vector<long double> wide(values.begin(),values.end());for(int k=0;k<dimension;k++)if(anchor[k]!=k)wide[anchor[k]]+=static_cast<long double>(volumes[variables[k].second])/volumes[variables[anchor[k]].second]*values[k];return std::vector<double>(wide.begin(),wide.end());};
    std::vector<saltchem::Inventory> state;std::vector<saltchem::AqueousState> seeds(n);std::vector<unsigned char> haveSeed(n);
    double exchangeResidual=0;
    auto evaluate=[&](const std::vector<double>& proposed,std::vector<saltchem::Inventory>& trial){
        std::vector<saltchem::MineralVector> extents(n);for(int k=0;k<dimension;k++)extents[variables[k].first][variables[k].second]=proposed[k];
        trial.resize(n);
        std::vector<double> porosity(n);
        for(int i=0;i<n;i++){
            trial[i]=saltchem::Reactions::changed(problems[i].initial,extents[i]);long double waterVolume=initialVolume[i];
            for(int m=0;m<3;m++){if(trial[i].mineral[m]<0)throw StepFailure("Local implicit reaction crossed mineral exhaustion");waterVolume+=static_cast<long double>(extents[i][m])*volumes[m];}
            if(!(waterVolume>1e-14)||(!retained[i]&&saltchem::aqueousMass(trial[i])<=0)){std::ostringstream message;message<<std::setprecision(17)<<"Local liquid volume floor: region="<<task.regions[i]<<" cell="<<g.waterRegions[task.regions[i]].parent<<" fraction="<<waterVolume<<" mass="<<saltchem::aqueousMass(trial[i]);throw StepFailure(message.str());}
            exchange.waterRegions[i].fraction=double(waterVolume);porosity[i]=double(waterVolume/problems[i].rates.geometry.bulkVolume);
        }
        auto reacted=trial;std::vector<double> density=initialDensity,targetMass(n);bool volumeConverged=false;
        for(int volumeIteration=0;volumeIteration<32;volumeIteration++){
            trial=reacted;for(int i=0;i<n;i++)targetMass[i]=density[i]*exchange.waterRegions[i].fraction;
            redistributeReactionWater(exchange,retained,targetMass,trial);
            // The active network is fixed during this nonlinear solve.
            // Crossing a representation threshold in a trial volume must
            // not discontinuously switch off its remaining water supply.
            auto transport=exchangeWater(exchange,trial,diffusivity,dt,0);stats[0].newtonIterations+=transport.iterations;exchangeResidual=transport.residual;
            double error=0;
            for(int i=0;i<n;i++)if(retained[i]){
                double water;auto properties=worker(g.waterRegions[task.regions[i]].parent).speciation.conserved(trial[i],water,haveSeed[i]?&seeds[i]:nullptr);
                stats[i].speciationIterations+=properties.iterations;seeds[i]=properties;haveSeed[i]=1;
                error=std::max(error,std::abs(properties.density/density[i]-1));density[i]=properties.density;
            }
            if(error<=std::max(1e-13,.1*options.relativeTolerance)){volumeConverged=true;break;}
        }
        if(!volumeConverged)throw StepFailure("Local water density-volume closure did not converge");
        std::vector<saltchem::MineralVector> rates(n);
        for(int i=0;i<n;i++){
            if(!(trial[i].waterEquivalent>0)||std::any_of(trial[i].aqueous.begin(),trial[i].aqueous.end(),[](double v){return v<0||!std::isfinite(v);}))throw StepFailure("Local implicit water composition crossed positivity");
            if(map.patches[task.regions[i]].empty())continue;
            auto potential=trial[i];
            // Evaluate the one-sided kinetic flux at exhaustion. The extent
            // residual enforces the available-stock bound below, so a finite
            // crystal can disappear within this implicit interval.
            for(int m=0;m<3;m++)if(problems[i].rates.area[m]>0&&potential.mineral[m]==0)potential.mineral[m]=1e-30;
            auto output=worker(g.waterRegions[task.regions[i]].parent).reactions.evaluate(potential,problems[i].rates,haveSeed[i]?&seeds[i]:nullptr,false,porosity[i]);
            rates[i]=output.rates;seeds[i]=output.aqueous;haveSeed[i]=1;stats[i].rateEvaluations++;stats[i].speciationIterations+=output.aqueous.iterations;
        }
        std::vector<double> residual(dimension);for(int k=0;k<dimension;k++){auto [i,m]=variables[k];residual[k]=proposed[k]-std::min(problems[i].initial.mineral[m],dt*rates[i][m]);}return residual;
    };
    auto extentBudget=[&](int k,const std::vector<double>& values,const std::vector<saltchem::Inventory>& trial){
        auto [i,m]=variables[k];const auto& state=trial[i];double absolute=options.absoluteTolerance*problems[i].initial.waterEquivalent;
        double budget=absolute+options.relativeTolerance*std::max(std::abs(values[k]),state.mineral[m]);
        for(int j=0;j<7;j++)if(saltchem::stoichiometry[m][j])budget=std::min(budget,(absolute+options.relativeTolerance*state.aqueous[j])/(3*saltchem::stoichiometry[m][j]));
        if(saltchem::hydration[m]>0)budget=std::min(budget,(saltchem::waterMolarMass*absolute+options.relativeTolerance*state.waterEquivalent)/(3*saltchem::hydration[m]*saltchem::waterMolarMass));
        return std::max(budget,1e-30);
    };
    auto norm=[&](const std::vector<double>& residual,const std::vector<double>& values,const std::vector<saltchem::Inventory>& trial){
        double maximum=0;for(int k=0;k<dimension;k++)maximum=std::max(maximum,std::abs(residual[k])/extentBudget(k,values,trial));return maximum;
    };
    auto stateError=[&](const std::vector<saltchem::Inventory>& a,const std::vector<saltchem::Inventory>& b){
        double maximum=0;for(int i=0;i<n;i++){double absolute=options.absoluteTolerance*problems[i].initial.waterEquivalent;
            for(int j=0;j<7;j++)maximum=std::max(maximum,std::abs(a[i].aqueous[j]-b[i].aqueous[j])/(absolute+options.relativeTolerance*std::max(std::abs(a[i].aqueous[j]),std::abs(b[i].aqueous[j]))));
            for(int m=0;m<3;m++)maximum=std::max(maximum,std::abs(a[i].mineral[m]-b[i].mineral[m])/(absolute+options.relativeTolerance*std::max(std::abs(a[i].mineral[m]),std::abs(b[i].mineral[m]))));
            maximum=std::max(maximum,std::abs(a[i].waterEquivalent-b[i].waterEquivalent)/(saltchem::waterMolarMass*absolute+options.relativeTolerance*std::max(a[i].waterEquivalent,b[i].waterEquivalent)));
        }return maximum;
    };
    const double requestedInterval=dt;double convergenceError=INFINITY;
    auto solveInterval=[&](double interval,const std::vector<double>& guess){
    dt=interval;x=guess;std::fill(haveSeed.begin(),haveSeed.end(),0);
    f=evaluate(x,state);bool ready=false,converged=false,preciseJacobian=false;convergenceError=INFINITY;
    for(int iteration=0;iteration<64;iteration++){
        double current=norm(f,x,state);if(current<=1){converged=true;convergenceError=current;break;}
        bool freshJacobian=false;
        if(!ready){
            for(int k=0;k<dimension;k++){
                auto [i,m]=variables[k];double water=std::max(problems[i].initial.waterEquivalent,state[i].waterEquivalent);
                double delta=preciseJacobian?std::max(.1*extentBudget(k,x,state),64*std::numeric_limits<double>::epsilon()*std::max(water,std::abs(x[k]))):1e-6*std::max(std::min(problems[i].initial.mineral[m],water),1e-6*water);bool done=false;
                for(int attempt=0;attempt<16&&!done;attempt++,delta*=.1)for(double sign:{1.,-1.}){
                    auto probe=x;probe[k]+=sign*delta;if(anchor[k]!=k)probe[anchor[k]]-=sign*delta*volumes[m]/volumes[variables[anchor[k]].second];std::vector<saltchem::Inventory> ignored;
                    try{auto values=evaluate(probe,ignored);for(int row=0;row<dimension;row++){double linear=row==k?1.:0.;if(anchor[k]!=k&&row==anchor[k])linear=-volumes[m]/volumes[variables[anchor[k]].second];jac[row*dimension+k]=linear+((values[row]-probe[row])-(f[row]-x[row]))/(sign*delta);}done=true;break;}catch(const std::exception&){}
                }
                if(!done){if(f[k]!=0)throw StepFailure("Local reaction-exchange derivative reached its admissible boundary");for(int row=0;row<dimension;row++){double linear=row==k?1.:0.;if(anchor[k]!=k&&row==anchor[k])linear=-volumes[m]/volumes[variables[anchor[k]].second];jac[row*dimension+k]=linear;}}
            }
            ready=true;
            freshJacobian=true;
        }
        auto rhs=f;for(double& v:rhs)v=-v;if(!solveDense(jac,rhs,direction)){ready=false;if(iteration%2==0)continue;throw StepFailure("Singular local reaction-exchange derivative");}
        direction=toPhysical(direction);
        const double correction=norm(direction,x,state);
        // Globalisation uses row-scaled kinetic equations. Otherwise one
        // almost-closed pore's amplified rate roundoff can block improvements
        // in every well-resolved neighbouring reservoir.
        auto merit=[&](const std::vector<double>& residual){auto scaled=residual;for(int k=0;k<dimension;k++)scaled[k]/=std::max(1.,std::abs(jac[k*dimension+k]));return norm(scaled,x,state);};
        const double currentMerit=merit(f);
        if(std::getenv("SALT_DUAL_NEWTON_TRACE")){std::cerr<<"newton "<<iteration<<" residual "<<current<<" correction "<<correction<<'\n';for(int k=0;k<dimension;k++){auto [i,m]=variables[k];std::cerr<<k<<" region="<<task.regions[i]<<" phase="<<m<<" water="<<problems[i].initial.waterEquivalent<<" volume="<<initialVolume[i]<<" stock="<<problems[i].initial.mineral[m]<<" x="<<x[k]<<" f="<<f[k]<<" d="<<direction[k]<<" diag="<<jac[k*dimension+k]<<'\n';}}
        // Stiff equilibria amplify a tiny state error into a large rate
        // residual. Verify a fresh-Jacobian correction and all inventories,
        // using the configured tolerances rather than shrinking physical time.
        if(iteration>0&&correction<=.5){
            if(!freshJacobian||!preciseJacobian){preciseJacobian=true;ready=false;continue;}
            auto proposed=x;std::vector<double> projected(dimension);for(int k=0;k<dimension;k++){auto [i,m]=variables[k];proposed[k]=std::min(problems[i].initial.mineral[m],x[k]+direction[k]);projected[k]=proposed[k]-x[k];}
            if(norm(projected,x,state)>=.1*correction){
                std::vector<saltchem::Inventory> trial;
                try{auto residual=evaluate(proposed,trial);double change=stateError(state,trial);if(change<=.5){x=std::move(proposed);f=std::move(residual);state=std::move(trial);convergenceError=std::max(correction,change);converged=true;break;}}
                catch(const std::exception&){}
            }
        }
        bool accepted=false;
        std::string trialFailure;
        double maximumAlpha=1;
        for(int k=0;k<dimension;k++){if(!std::isfinite(direction[k]))throw StepFailure("Nonfinite local Newton step");if(direction[k]!=0){int i=variables[k].first;maximumAlpha=std::min(maximumAlpha,10*std::max({problems[i].initial.waterEquivalent,state[i].waterEquivalent,std::abs(x[k])})/std::abs(direction[k]));}}
        for(int backtrack=0;backtrack<32;backtrack++){
            double alpha=std::ldexp(maximumAlpha,-backtrack);
            auto proposed=x;for(int k=0;k<dimension;k++){auto [i,m]=variables[k];proposed[k]=std::min(problems[i].initial.mineral[m],x[k]+alpha*direction[k]);}std::vector<saltchem::Inventory> trial;
            try{
                auto residual=evaluate(proposed,trial);double error=merit(residual);
                if(std::getenv("SALT_DUAL_NEWTON_TRACE"))std::cerr<<" alpha="<<alpha<<" error="<<error<<'\n';
                if(error<currentMerit*(1-1e-4*alpha)||norm(residual,proposed,trial)<=1){
                    std::vector<double> difference(dimension);for(int k=0;k<dimension;k++)difference[k]=proposed[k]-x[k];difference=toBasis(difference);double denominator=0;for(double value:difference)denominator+=value*value;
                    if(denominator>1e-60)for(int row=0;row<dimension;row++){double update=residual[row]-f[row];for(int k=0;k<dimension;k++)update-=jac[row*dimension+k]*difference[k];for(int k=0;k<dimension;k++)jac[row*dimension+k]+=update*difference[k]/denominator;}
                    x=std::move(proposed);f=std::move(residual);state=std::move(trial);ready=alpha>.125&&iteration%8!=7;preciseJacobian=false;accepted=true;break;
                }
            }catch(const std::exception& error){trialFailure=error.what();}
        }
        if(!accepted){
            // At a phase-appearance/exhaustion corner the Newton direction
            // can point outside the feasible cone. A damped kinetic residual
            // step preserves the stock bound while rebuilding that Jacobian.
            std::vector<double> fallback(dimension);for(int k=0;k<dimension;k++)fallback[k]=-f[k]/std::max(1.,std::abs(jac[k*dimension+k]));fallback=toPhysical(fallback);
            double maximumFallback=1;for(int k=0;k<dimension;k++){double d=fallback[k];if(d!=0){int i=variables[k].first;maximumFallback=std::min(maximumFallback,10*std::max({problems[i].initial.waterEquivalent,state[i].waterEquivalent,std::abs(x[k])})/std::abs(d));}}
            for(int backtrack=0;backtrack<32;backtrack++){
                double alpha=std::ldexp(maximumFallback,-backtrack);
                auto proposed=x;for(int k=0;k<dimension;k++){auto [i,m]=variables[k];proposed[k]=std::min(problems[i].initial.mineral[m],x[k]+alpha*fallback[k]);}
                std::vector<saltchem::Inventory> trial;
                try{auto residual=evaluate(proposed,trial);double error=merit(residual);if(error<currentMerit*(1-1e-4*alpha)||norm(residual,proposed,trial)<=1){x=std::move(proposed);f=std::move(residual);state=std::move(trial);ready=false;preciseJacobian=false;accepted=true;break;}}
                catch(const std::exception& error){trialFailure=error.what();}
            }
            if(accepted)continue;
            if(ready&&iteration%2==0){ready=false;continue;}
            std::ostringstream message;message<<"Local reaction-exchange solve requested a shorter interval: iteration="<<iteration<<" residual="<<current<<" detail="<<trialFailure;throw StepFailure(message.str());
        }
    }
    if(!converged){convergenceError=norm(f,x,state);if(convergenceError>1)throw StepFailure("Local reaction-exchange state error exceeded tolerance");}
    };
    try{solveInterval(requestedInterval,std::vector<double>(dimension));}
    catch(const StepFailure& first){
        // Continuation supplies Newton guesses for the same initial inventory.
        // Only the root for the original requested interval is committed.
        std::vector<double> previous(dimension);double solved=0,trial=.25*requestedInterval;bool reached=false;std::string detail=first.what();
        for(int stage=0;stage<40;stage++){
            if(trial-solved<std::max(1e-15,requestedInterval*1e-10))break;
            try{solveInterval(trial,previous);previous=x;solved=trial;if(trial==requestedInterval){reached=true;break;}trial=std::min(requestedInterval,2*trial);}
            catch(const StepFailure& error){detail=error.what();trial=solved>0?.5*(solved+trial):.25*trial;}
        }
        if(!reached)throw StepFailure("Local reaction-exchange continuation requested a shorter interval: "+detail);
    }
    const double finalNorm=norm(f,x,state);std::vector<saltchem::MineralVector> extents(n);for(int k=0;k<dimension;k++)extents[variables[k].first][variables[k].second]=x[k];
    for(int i=0;i<n;i++){
        stats[i].accepted++;
        for(int m=0;m<3;m++){double extent=extents[i][m];stats[i].dissolved[m]=std::max(0.,extent);stats[i].precipitated[m]=std::max(0.,-extent);}
        result.push_back(regionResult(task.regions[i],g,map,cells,volumes,problems[i].initial,state[i],stats[i],&extents[i]));
    }
    result.front().exchangeResidual=exchangeResidual;result.front().nonlinearResidualRatio=finalNorm;
    result.front().nonlinearStateErrorRatio=convergenceError;
    return result;
}
}
