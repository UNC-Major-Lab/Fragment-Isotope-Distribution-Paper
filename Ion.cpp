//
// Created by Mike Lafferty on 9/16/16.
//

#include <ostream>
#include <iomanip>
#include "Ion.h"

Ion::Ion(OpenMS::AASequence seq, OpenMS::Residue::ResidueType type, OpenMS::Int charge) {
    this->sequence = seq;
    this->type = type;
    this->charge = charge;
    this->formula = seq.getFormula(type, charge);
    this->monoWeight = seq.getMonoWeight(type, charge);
    this->monoMz = this->monoWeight / charge;
}

std::vector<Ion> Ion::generateFragmentIons(double minMz, double maxMz) {
    std::vector<Ion> ionList;
    //generate b-ions
    //loop through each b-ion
    for (int i = 1; i <= sequence.size(); ++i) {
        //loop through each charge
        for (int z = 1; z < charge; ++z) {
            //add b-ion to list
            double mz = sequence.getPrefix(i).getMonoWeight(OpenMS::Residue::BIon, z) / z;
            if (mz >= minMz && mz <= maxMz) {
                ionList.push_back(Ion(sequence.getPrefix(i), OpenMS::Residue::BIon, z));
            }
        }
    }
    //generate y-ions
    //loop through each y-ion
    for (int i = 1; i <= sequence.size() - 1; ++i) {
        //loop through each charge
        for (int z = 1; z < charge; ++z) {
            //add y-ion to list
            double mz = sequence.getSuffix(i).getMonoWeight(OpenMS::Residue::YIon, z) / z;
            if (mz >= minMz && mz <= maxMz) {
                ionList.push_back(Ion(sequence.getSuffix(i), OpenMS::Residue::YIon, z));
            }
        }
    }
    //generate precursor ion
    double mz = sequence.getMonoWeight(OpenMS::Residue::Full, charge) / charge;
    if (mz >= minMz && mz <= maxMz) {
        ionList.push_back(Ion(sequence, OpenMS::Residue::Full, charge));
    }

    return ionList;
}

std::string Ion::getIonType() {
    std::string ion_type = (type == OpenMS::Residue::ResidueType::BIon) ? "B" : "Y";
    ion_type += std::to_string(sequence.size());
    for (int i = 0; i < charge; ++i) ion_type += "+";
    return ion_type;
}

std::string Ion::getIonName() {
    return getIonType() + " " + sequence.toUnmodifiedString();
}

std::ostream& operator<<(std::ostream &strm, const Ion &ion) {
    //print all ion information
    strm << "Ion";
    strm << " | type: " << std::setw(2) << ion.type;
    strm << " | charge: " << std::setw(2) << ion.charge;
    strm << " | weight: " << std::setw(8) << ion.monoWeight;
    strm << " | mz: " << std::setw(8) << ion.monoWeight / ion.charge;
    strm << " | seq: " << ion.sequence;
    //return output stream
    return strm;
}
