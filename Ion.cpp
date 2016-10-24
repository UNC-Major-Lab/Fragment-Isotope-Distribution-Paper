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
    this->mz = this->monoWeight / charge;
}

void Ion::generateFragmentIons(std::vector<Ion> &ionList, const OpenMS::AASequence pepSeq, const OpenMS::Int pepCharge) {
    //generate b-ions
    //loop through each b-ion
    for (int i = 1; i <= pepSeq.size() - 1; ++i) {
        //loop through each charge
        for (int charge = 1; charge <= pepCharge; ++charge) {
            //add b-ion to list
            ionList.push_back(Ion(pepSeq.getPrefix(i), OpenMS::Residue::BIon, charge));
        }
    }
    //generate y-ions
    //loop through each y-ion
    for (int i = 1; i <= pepSeq.size() - 1; ++i) {
        //loop through each charge
        for (int charge = 1; charge <= pepCharge; ++charge) {
            //add b-ion to list
            ionList.push_back(Ion(pepSeq.getSuffix(i), OpenMS::Residue::YIon, charge));
        }
    }
    //generate precursor ion
    ionList.push_back(Ion(pepSeq, OpenMS::Residue::Full, pepCharge));
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
