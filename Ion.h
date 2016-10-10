// Ion class to hold fragment ions from precursor peptides.
//
// Created by Mike Lafferty on 9/16/16.
//

#ifndef SPECOPS_ION_H
#define SPECOPS_ION_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

class Ion {
public:
    OpenMS::Residue::ResidueType type;
    OpenMS::AASequence sequence;
    OpenMS::Int charge;
    OpenMS::EmpiricalFormula formula;
    double monoWeight;

    /**
     * Ion constructor. Includes public members for sequence, type, charge, molecular formula,
     * and monoisotopic weight (charge included) in Daltons.
     * @param seq the amino acid sequence of the ion
     * @param type the type of ion (y-ion,b-ion,full,...) See OpenMS::Residue
     * @param charge the ion charge
     * @return a constructed Ion object
     */
    Ion(OpenMS::AASequence seq, OpenMS::Residue::ResidueType type, OpenMS::Int charge);

    /**
     * Generates a list of fragment y- and b-ions from an input peptide sequence and charge and addes them
     * to a vector of Ion objects. Ions will be created of charge 1 and all charges up to and including the
     * precursor peptide charge. The precursor peptide is also added to the list of fragment ions.
     * @param ionList a vector of Ion objects. If not empty, generated ions will be appended to the end.
     * @param pepSeq the amino acid sequence of the precursor peptide to generate fragment ions from
     * @param pepCharge the charge of the precursor peptide
     */
    static void generateFragmentIons(std::vector<Ion> &ionList, const OpenMS::AASequence pepSeq,
                                     const OpenMS::Int pepCharge);

    /**
     * Function to send an ion to standard output stream. All ion members reported in addition
     * to the ion mz.
     * @param strm output stream
     * @param ion ion to report
     * @return output stream with ion information
     */
    friend std::ostream& operator<<(std::ostream &strm, const Ion &ion);
};

#endif //SPECOPS_ION_H
