// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include <array>
// Prescribed spatially uniform conditions; rows are (time_s temperature_C
// absolute_pressure_atm). Linear interpolation, constant outside the table.
class ThermodynamicSchedule {
    Foam::List<Foam::vector> table;
    double temperature,pressure;
public:
    explicit ThermodynamicSchedule(const Foam::dictionary& cfg):
        temperature(cfg.lookupOrDefault<Foam::scalar>("temperatureC",25)),
        pressure(cfg.lookupOrDefault<Foam::scalar>("absolutePressureAtm",1)) {
        if(cfg.found("thermodynamicSchedule"))cfg.lookup("thermodynamicSchedule")>>table;
        auto validate=[](double t,double p){if(!std::isfinite(t)||t+273.15<=228||!std::isfinite(p)||p<=0)throw std::invalid_argument("Invalid temperatureC or absolutePressureAtm");};
        validate(temperature,pressure);
        for(Foam::label i=0;i<table.size();i++){
            validate(table[i].y(),table[i].z());
            if(!std::isfinite(table[i].x())||(i&&table[i].x()<=table[i-1].x()))throw std::invalid_argument("Thermodynamic schedule times must increase strictly");
        }
    }
    std::array<double,2> at(double time)const {
        if(table.empty())return {{temperature,pressure}};
        if(time<=table.first().x())return {{table.first().y(),table.first().z()}};
        for(Foam::label i=1;i<table.size();i++)if(time<=table[i].x()){
            const double w=(time-table[i-1].x())/(table[i].x()-table[i-1].x());
            return {{table[i-1].y()+w*(table[i].y()-table[i-1].y()),table[i-1].z()+w*(table[i].z()-table[i-1].z())}};
        }return {{table.last().y(),table.last().z()}};
    }
    double nextKnot(double time,double end)const {
        for(const auto& row:table)if(row.x()>time+1e-12)end=std::min(end,double(row.x()));
        return end;
    }
};
