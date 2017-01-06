//
// Created by Dennis Goldfarb on 1/6/17.
//

#include <iostream>
#include <map>

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

#include "FASTAParser.h"

using namespace OpenMS;

static const OpenMS::ElementDB* elementDB = OpenMS::ElementDB::getInstance();

void usage()
{
    std::cout << "CreateAveragineModel fasta_path min_peptide_length max_peptide_length" << std::endl;
}

int main(int argc, const char ** argv) {

    if (argc != 4)
    {
        usage();
    }

    std::string fasta_path = argv[1];
    int min_peptide_length = atoi(argv[2]);
    int max_peptide_length = atoi(argv[3]);

    EmpiricalFormula ef;
    double count = 0;

    FASTAParser parser(fasta_path, 1e10, min_peptide_length, max_peptide_length);

    for (auto itr = parser.begin(); itr != parser.end(); ++itr)
    {
        ef += itr->getFormula();
        count+= itr->size();
    }

    for (auto itr = ef.begin(); itr != ef.end(); ++itr)
    {
        std::cout << itr->first->getSymbol() << "\t" << itr->second/count << std::endl;
    }

    return 0;
}