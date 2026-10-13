// SPDX-License-Identifier: GPL-3.0-or-later
#include "BrineReactions.hpp"
#include <fstream>
#include <iostream>
using namespace saltchem;
int main(int argc,char** argv){try{
    const std::string db=argc>1?argv[1]:"src/chemistry/pitzer.dat";
    std::ifstream in(db);if(!in)throw std::runtime_error("Missing database");std::string text((std::istreambuf_iterator<char>(in)),{});
    std::replace(text.begin(),text.end(),';','\n');DatabaseContract::validate(text);int rejected=0;
    auto reject=[&](std::string data){try{DatabaseContract::validate(data);}catch(const std::exception&){rejected++;return;}throw std::runtime_error("Incompatible database was accepted");};
    auto replacement=[&](const std::string& a,const std::string& b){std::string s=text;const auto at=s.find(a);if(at==std::string::npos)throw std::runtime_error("Missing mutation target: "+a);s.replace(at,a.size(),b);return s;};
    reject(replacement("KMgCl3:6H2O = K+ + Mg+2 + 3 Cl- + 6 H2O","KMgCl3:6H2O = K+ + Mg+2 + 2 Cl- + 6 H2O"));
    reject(replacement("CO3-2 + H+ = HCO3-","CO3-2 + 2 H+ = HCO3-"));
    reject(replacement("SOLUTION_SPECIES","SOLUTION_SPECIES\nNa+ + OH- = NaOH\nlog_k 0"));
    reject(replacement("PITZER","PITZER\n-macinnes false"));
    reject(replacement("PITZER","PITZER\n-use_etheta false"));
    reject(replacement("PITZER","PITZER\n-redox true"));
    PitzerModel model(db,650,1);auto w=model.waterState(.8);if(!w.temperatureCapped||!w.pressureRaised||w.waterTemperatureK!=623.15)throw std::runtime_error("Missing water diagnostics");
    bool strict=false;try{model.setStrictWater(true);}catch(const std::domain_error&){strict=true;}if(!strict)throw std::runtime_error("Strict water cap policy not enforced");
    model.setConditions(423.15,1);model.setStrictWater(true);strict=false;try{Speciation(model).solve(Totals{{1,0,0,0,1,0,0}});}catch(const std::domain_error&){strict=true;}if(!strict)throw std::runtime_error("Strict effective pressure policy not enforced");
    RateParameters p;p.geometry.enabled=true;p.geometry.closureWidth=0;bool bad=false;try{Reactions::validateParameters(p);}catch(const std::invalid_argument&){bad=true;}if(!bad)throw std::runtime_error("Invalid closure width accepted");
    std::cout<<"database_rejections="<<rejected<<" water_diagnostics=passed strict_policy=passed closure_config=passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
