// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrineSystem.hpp"
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <cstdlib>
namespace saltchem {
class DatabaseContract {
    using Side=std::map<std::string,double>;
    static Side side(const std::string& text){
        std::istringstream in(text);std::string token;double multiplier=1;Side out;
        while(in>>token){
            if(token=="+")continue;
            char* end=nullptr;const double number=std::strtod(token.c_str(),&end);
            if(end!=token.c_str()&&*end=='\0'){
                if(!std::isfinite(number)||number<=0)throw std::runtime_error("Invalid reaction coefficient");
                multiplier=number;continue;
            }
            out[token]+=multiplier;multiplier=1;
        }
        return out;
    }
    static std::pair<Side,Side> sides(const std::string& line){
        const auto at=line.find('=');
        if(at==std::string::npos)throw std::runtime_error("Missing reaction equation");
        return {side(line.substr(0,at)),side(line.substr(at+1))};
    }
    static Side net(const std::pair<Side,Side>& reaction){
        Side n=reaction.second;for(const auto& v:reaction.first)n[v.first]-=v.second;
        for(auto p=n.begin();p!=n.end();)if(p->second==0)p=n.erase(p);else ++p;
        return n;
    }
    static bool known(const Side& s){for(const auto& v:s)if(v.first!="H2O"&&speciesIndex(v.first)<0)return false;return true;}
    static const std::array<const char*,6>& aqueous(){
        static const std::array<const char*,6> a{{"H2O = OH- + H+","CO3-2 + H+ = HCO3-","CO3-2 + 2 H+ = CO2 + H2O","SO4-2 + H+ = HSO4-","Mg+2 + H2O = MgOH+ + H+","Mg+2 + CO3-2 = MgCO3"}};return a;
    }
public:
    static int aqueousIndex(const std::string& line){
        const auto reaction=sides(line);const auto n=net(reaction);
        if(n.empty()&&reaction.first.size()==1)return speciesIndex(reaction.first.begin()->first);
        const int index[6]={OH,Bicarbonate,CO2,Bisulfate,MgOH,MgCO3};
        for(int i=0;i<6;i++)if(n==net(sides(aqueous()[i])))return index[i];
        if(known(reaction.first)||known(reaction.second))throw std::runtime_error("Reaction outside the compiled aqueous network: "+line);
        return -1;
    }
    static void phase(int mineral,const std::string& line){
        static const char* formula[7]={"NaCl","KCl","KMgCl3:6H2O","K2MgCa2(SO4)4:2H2O","CaSO4:2H2O","MgCl2:6H2O","CaCO3"};
        const auto reaction=sides(line);Side expected;
        for(int i=0;i<nComponents;i++)if(stoichiometry[mineral][i])expected[speciesNames[i]]=stoichiometry[mineral][i];
        if(hydration[mineral])expected["H2O"]=hydration[mineral];
        if(reaction.first!=Side{{formula[mineral],1}}||reaction.second!=expected)
            throw std::runtime_error("Mineral stoichiometry does not match compiled network: "+std::string(mineralNames[mineral]));
    }
    static void validate(const std::string& database){
        std::istringstream lines(database);std::string line,section;
        std::array<bool,nSpecies> speciesSeen{},volumeSeen{};
        std::array<bool,nMinerals> phaseSeen{};
        const std::map<std::string,std::pair<std::string,double>> masters={
            {"Na",{"Na+",22.9898}},{"K",{"K+",39.0983}},{"Ca",{"Ca+2",40.08}},{"Mg",{"Mg+2",24.305}},
            {"Cl",{"Cl-",35.453}},{"S",{"SO4-2",32.064}},{"C",{"CO3-2",12.0111}},{"H",{"H+",1.008}},{"O",{"H2O",16.00}}};
        std::set<std::string> masterSeen;
        int species=-1,mineral=-1;int number=0;
        const std::set<std::string> sections={"SOLUTION_MASTER_SPECIES","SOLUTION_SPECIES","PHASES","PITZER","GAS_BINARY_PARAMETERS","EXCHANGE_MASTER_SPECIES","EXCHANGE_SPECIES","SURFACE_MASTER_SPECIES","SURFACE_SPECIES","RATES","MEAN_GAMMAS","END"};
        try{while(std::getline(lines,line)){
            ++number;line=line.substr(0,line.find('#'));std::istringstream in(line);std::string first;if(!(in>>first))continue;
            if(sections.count(first)){section=first;species=-1;mineral=-1;continue;}
            if(section=="SOLUTION_MASTER_SPECIES"){
                const auto match=masters.find(first);
                if(match!=masters.end()){
                    std::string speciesName,formula;double alkalinity,weight;
                    if(!(in>>speciesName>>alkalinity>>formula>>weight)||speciesName!=match->second.first||!std::isfinite(weight)||std::abs(weight-match->second.second)>1e-10)
                        throw std::runtime_error("Master species or atomic mass differs from compiled inventory basis: "+first);
                    masterSeen.insert(first);
                }
            }else if(section=="SOLUTION_SPECIES"){
                if(line.find('=')!=std::string::npos){species=aqueousIndex(line);if(species>=0)speciesSeen[species]=true;continue;}
                if(species>=0){
                    if(first=="-Vm")volumeSeen[species]=true;
                    else if(first!="log_k"&&first!="-log_k"&&first!="delta_h"&&first!="-delta_h"&&first.rfind("-analytic",0)!=0&&first!="-dw"&&first!="-viscosity")
                        throw std::runtime_error("Unsupported active aqueous option: "+first);
                }
            }else if(section=="PHASES"){
                if(!std::isspace(static_cast<unsigned char>(line[0]))){mineral=mineralIndex(first);continue;}
                if(mineral>=0){
                    if(line.find('=')!=std::string::npos){phase(mineral,line);phaseSeen[mineral]=true;}
                    else if(first!="Vm"&&first!="-Vm"&&first!="log_k"&&first!="-log_k"&&first!="delta_h"&&first!="-delta_h"&&first.rfind("-analytic",0)!=0)
                        throw std::runtime_error("Unsupported active phase option: "+first);
                }
            }else if(section=="PITZER"&&(first=="-macinnes"||first=="-use_etheta"||first=="-redox")){
                std::string value;in>>value;const bool enabled=value.empty()||value=="true"||value=="t";
                if(enabled!=(first!="-redox"))throw std::runtime_error("Pitzer switch differs from compiled convention: "+first);
            }
        }
        for(int i=0;i<nSpecies;i++){
            if(!speciesSeen[i])throw std::runtime_error("Missing aqueous species: "+std::string(speciesNames[i]));
            // H+ is PHREEQC's reference volume; pitzer.dat explicitly uses the
            // zero default for MgOH+. Other active species require volume data.
            if(i!=H&&i!=MgOH&&!volumeSeen[i])throw std::runtime_error("Missing aqueous volume parameters: "+std::string(speciesNames[i]));
        }
        if(masterSeen.size()!=masters.size())throw std::runtime_error("Missing master species or atomic masses");
        for(int i=0;i<nMinerals;i++)if(!phaseSeen[i])throw std::runtime_error("Missing phase reaction: "+std::string(mineralNames[i]));
        }catch(const std::exception& e){throw std::runtime_error("Pitzer database contract near line "+std::to_string(number)+": "+e.what());}
    }
};
}
