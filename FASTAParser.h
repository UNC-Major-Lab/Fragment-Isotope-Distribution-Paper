//
// Created by Dennis Goldfarb on 1/6/17.
//

#ifndef FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_FASTAPARSER_H
#define FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_FASTAPARSER_H

#include <iostream>
#include <set>

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

using namespace OpenMS;

class FASTAParser {

public:

    FASTAParser(std::string fasta_path, float max_mass, UInt min_peptide_length, UInt max_peptide_length) :
            fasta_path(fasta_path), MAX_MASS(max_mass), MIN_PEPTIDE_LENGTH(min_peptide_length), MAX_PEPTIDE_LENGTH(max_peptide_length)
    {
        digestFASTA();
    }

    typedef std::set<AASequence>::iterator iterator;
    typedef std::set<AASequence>::iterator Iterator;

    inline Iterator begin() { return unique_peptides.begin(); }
    inline Iterator end()   { return unique_peptides.end(); }

private:
    bool isValidPeptide(AASequence& pep);
    void digestFASTA();
    void digestProtein(FASTAFile::FASTAEntry& protein, EnzymaticDigestion& digestor);

    std::string fasta_path;
    std::set<AASequence> unique_peptides;

    UInt MIN_PEPTIDE_LENGTH;
    UInt MAX_PEPTIDE_LENGTH;
    float MAX_MASS:
};


#endif //FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_FASTAPARSER_H
