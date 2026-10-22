// SPDX-License-Identifier: GPL-3.0-or-later
#include "DualMineralGeometry.hpp"
#include <iostream>
int main(int argc,char** argv){
    if(argc!=3){std::cerr<<"Usage: makeDualGeometry output.bin bare|seeded|mixed\n";return 2;}
    saltdual::Geometry g(64,40,4,5e-6,20e-6,1e-5);
    if(std::string(argv[2])=="mixed"){
        g.circle(0,160e-6,65e-6,20e-6);
        g.circle(1,220e-6,115e-6,28e-6);
        g.circle(2,100e-6,115e-6,28e-6);
    }else g.circle(2,160e-6,100e-6,40e-6);
    g.rebuild();int seeds=0;
    if(std::string(argv[2])=="seeded")for(int k=0;k<24;k++){
        double angle=2*saltgeom::pi*(k+.37)/24;
        g.circle(k%2,160e-6+40.5e-6*std::cos(angle),100e-6+40.5e-6*std::sin(angle),1.25e-6,64);
        seeds++;
    }else if(std::string(argv[2])!="bare"&&std::string(argv[2])!="mixed")return 2;
    g.rebuild();
    if(g.globalLiquidComponents!=1)throw std::runtime_error("Initial brine must remain globally connected");
    g.write(argv[1]);std::ofstream f(std::string(argv[1])+".coarse",std::ios::binary);
    for(auto a:g.fractions)f.write(reinterpret_cast<const char*>(a.data()),sizeof(a));
    std::cout<<"Main cells "<<g.nx*g.ny<<" geometry cells "<<g.cells.size()<<" wet segments "<<g.wet.size()<<" initial_seeds "<<seeds<<'\n';
}
