//
// Created by Dennis Goldfarb on 11/21/16.
//

#include <iostream>

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

using namespace OpenMS;

static const OpenMS::ElementDB* elementDB = OpenMS::ElementDB::getInstance();

std::set<AASequence> uniquePeptides;
std::map<int, int> sulfurs2count;

bool isValidPeptide(AASequence& pep) {
    String p = pep.toString();
    if (p.hasSubstring("U") || p.hasSubstring("B") || p.hasSubstring("Z") || p.hasSubstring("J") || p.hasSubstring("X"))
    {
        return false;
    }
    return true;
}

void outputDistribution()
{
    for (auto itr : sulfurs2count)
    {
        std::cout << itr.first << "\t" << itr.second << std::endl;
    }
}

void countSulfurs(AASequence& pep)
{
    int pep_s = pep.getFormula().getNumberOf(ElementDB::getInstance()->getElement("Sulfur"));

    if (sulfurs2count.find(pep_s) == sulfurs2count.end()) {
        sulfurs2count[pep_s] = 0;
    }

    sulfurs2count[pep_s]++;
}

void digestProtein(FASTAFile::FASTAEntry& protein, EnzymaticDigestion& digestor)
{
    static Size MIN_PEPTIDE_LENGTH = 5;
    static Size MAX_PEPTIDE_LENGTH = 100;

    std::vector<AASequence> peptides;
    digestor.digest(AASequence::fromString(protein.sequence), peptides);
    for (Size j = 0; j < peptides.size(); ++j)
    {
        if (peptides[j].size() >= MIN_PEPTIDE_LENGTH && peptides[j].size() <= MAX_PEPTIDE_LENGTH
            && isValidPeptide(peptides[j]) && uniquePeptides.find(peptides[j]) == uniquePeptides.end())
        {
            uniquePeptides.insert(peptides[j]);
            countSulfurs(peptides[j]);
        }
    }
}

void digestFASTA(std::string fasta_path)
{
    std::vector<FASTAFile::FASTAEntry> proteins;
    FASTAFile().load(fasta_path, proteins);

    EnzymaticDigestion digestor; // default parameters are fully tryptic with 0 missed cleavages

    for (Size i = 0; i < proteins.size(); ++i)
    {
        digestProtein(proteins[i], digestor);
    }
}

void usage()
{
    std::cout << "GetSulfurDistribution fasta_path" << std::endl;
}

int main(int argc, const char ** argv) {

    if (argc != 2)
    {
        usage();
    }

    digestFASTA(argv[1]);
    outputDistribution();

    return 0;
}