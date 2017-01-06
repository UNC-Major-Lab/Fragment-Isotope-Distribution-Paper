//
// Created by Dennis Goldfarb on 1/6/17.
//

#include "FASTAParser.h"

bool FASTAParser::isValidPeptide(AASequence& pep) {
    String p = pep.toString();
    if (p.hasSubstring("U") || p.hasSubstring("B") || p.hasSubstring("Z") || p.hasSubstring("J") || p.hasSubstring("X"))
    {
        return false;
    }
    return true;
}

void FASTAParser::digestProtein(FASTAFile::FASTAEntry& protein, EnzymaticDigestion& digestor)
{
    std::vector<AASequence> peptides;
    digestor.digest(AASequence::fromString(protein.sequence), peptides);
    for (Size j = 0; j < peptides.size(); ++j)
    {
        if (peptides[j].size() >= MIN_PEPTIDE_LENGTH && peptides[j].size() <= MAX_PEPTIDE_LENGTH
            && isValidPeptide(peptides[j]) && unique_peptides.find(peptides[j]) == unique_peptides.end())
        {
            unique_peptides.insert(peptides[j]);
        }
    }
}

void FASTAParser::digestFASTA()
{
    std::vector<FASTAFile::FASTAEntry> proteins;
    FASTAFile().load(fasta_path, proteins);

    EnzymaticDigestion digestor; // default parameters are fully tryptic with 0 missed cleavages

    for (Size i = 0; i < proteins.size(); ++i)
    {
        digestProtein(proteins[i], digestor);
    }
}