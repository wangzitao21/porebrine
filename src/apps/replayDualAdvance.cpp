// SPDX-License-Identifier: GPL-3.0-or-later
#include "DualMineralGeometry.hpp"
#include <iostream>
int main(int argc,char** argv){try{
    if(argc!=3)throw std::invalid_argument("Supply geometry and reaction-budget binary files");
    auto g=saltdual::Geometry::read(argv[1]);auto before=g.fractions;std::ifstream f(argv[2],std::ios::binary);auto target=before;
    for(auto& a:target)f.read(reinterpret_cast<char*>(a.data()),sizeof(a));std::uint32_t count=0;f.read(reinterpret_cast<char*>(&count),4);if(count>1000000)throw std::runtime_error("Invalid budget count");
    std::vector<saltdual::SurfaceBudget> budgets(count);
    for(auto& b:budgets){std::int32_t values[4];f.read(reinterpret_cast<char*>(values),sizeof(values));b.parent=values[0];b.mineral=values[1];b.liquidRegion=values[2];if(values[3]<0||values[3]>1000000)throw std::runtime_error("Invalid surface count");b.surfaces.resize(values[3]);f.read(reinterpret_cast<char*>(&b.change),8);f.read(reinterpret_cast<char*>(b.surfaces.data()),4*b.surfaces.size());}
    if(!f)throw std::runtime_error("Incomplete reaction budgets");bool passed=true;
    if(std::getenv("SALT_DUAL_GEOM_TRACE")){for(int m=0;m<3;m++)std::cerr<<"debt "<<m<<' '<<g.chemicalVolumeDebt[m]<<'\n';for(const auto& b:budgets)if(b.liquidRegion>=0&&!g.waterRegions[b.liquidRegion].bulk)std::cerr<<"retained budget region="<<b.liquidRegion<<" phase="<<b.mineral<<" change="<<b.change<<" water="<<g.waterRegions[b.liquidRegion].fraction<<" faces="<<b.surfaces.size()<<'\n';}
    try{g.advance(target,false,.2,&budgets);}catch(const std::exception& e){std::cerr<<e.what()<<'\n';passed=false;}
    std::cout<<std::setprecision(17);for(int m=0;m<3;m++){long double wanted=0,actual=0;for(std::size_t p=0;p<target.size();p++){wanted+=static_cast<long double>(target[p][m])-before[p][m];actual+=static_cast<long double>(g.fractions[p][m])-before[p][m];}std::cout<<"phase="<<m<<" wanted="<<double(wanted)<<" actual="<<double(actual)<<" error="<<double(actual-wanted)<<'\n';}
    return passed?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
