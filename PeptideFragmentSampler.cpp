//
// Created by Dennis Goldfarb on 9/28/16.
//

#include <iostream>
#include <fstream>
#include <random>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>


static const OpenMS::ResidueDB* residueDB = OpenMS::ResidueDB::getInstance();
static const OpenMS::ElementDB* elementDB = OpenMS::ElementDB::getInstance();

static std::string AMINO_ACIDS = "ADEFGHIKLNPQRSTVWY";
static std::string AMINO_ACIDS_SULFUR = "CM";
static std::string AMINO_ACIDS_SELENIUM = "U";

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dis_AA(0, AMINO_ACIDS.length()-1);
std::uniform_int_distribution<> dis_S(0, AMINO_ACIDS_SULFUR.length()-1);

const int max_depth = 6;

OpenMS::AASequence create_random_peptide_sequence(int peptide_length, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium) {
    OpenMS::AASequence random_peptide;

    // for insertion of sulfur containing amino acids in fragment
    for (int i = 0; i < num_sulfurs; ++i) {
        random_peptide += residueDB->getResidue(AMINO_ACIDS_SULFUR[dis_S(gen)]);
    }

    // for insertion of selenocysteines in fragment
    for (int i = 0; i < num_selenium; ++i) {
        random_peptide += residueDB->getResidue(AMINO_ACIDS_SELENIUM[0]);
    }

    // random amino acid insertion (non Sulfur and Selenium amino acids)
    for (int aa_index = 0; aa_index < peptide_length; ++aa_index) {
        random_peptide += residueDB->getResidue(AMINO_ACIDS[dis_AA(gen)]);
    }

    // for insertion of sulfur containing amino acids in fragment
    for (int i = 0; i < num_c_sulfurs; ++i) {
        random_peptide += residueDB->getResidue(AMINO_ACIDS_SULFUR[dis_S(gen)]);
    }

    // for insertion of selenocysteines in fragment
    for (int i = 0; i < num_c_selenium; ++i) {
        random_peptide += residueDB->getResidue(AMINO_ACIDS_SELENIUM[0]);
    }


    return random_peptide;
}

void create_fragments(OpenMS::AASequence &p, std::ofstream outfiles[max_depth][max_depth], int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium) {
    int num_fragments = p.size()-1;

    int tot_left_SSe = std::max(num_sulfurs + num_selenium,1);
    int tot_right_SSe = std::max(num_c_sulfurs + num_c_selenium,1);

    OpenMS::EmpiricalFormula precursor_ef = p.getFormula();

    for (int index = tot_left_SSe; index < num_fragments-tot_right_SSe; ++index)
    {
        OpenMS::EmpiricalFormula b_ion = p.getPrefix(index).getFormula(OpenMS::Residue::ResidueType::BIon);
        OpenMS::EmpiricalFormula y_ion = p.getPrefix(index).getFormula(OpenMS::Residue::ResidueType::YIon);

        for (int precursor_isotope = 0; precursor_isotope < max_depth; ++precursor_isotope)
        {
            std::vector<OpenMS::UInt> isolated_isotopes;
            isolated_isotopes.push_back(precursor_isotope);

            OpenMS::IsotopeDistribution b_id = b_ion.getConditionalFragmentIsotopeDist(precursor_ef, isolated_isotopes);
            OpenMS::IsotopeDistribution y_id = y_ion.getConditionalFragmentIsotopeDist(precursor_ef, isolated_isotopes);

            for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope && fragment_isotope < b_id.size(); ++fragment_isotope)
            {
                outfiles[precursor_isotope][fragment_isotope] << b_id.getContainer()[fragment_isotope].second
                                                              << "\t" << b_ion.getMonoWeight()
                                                              << "\t" << p.getMonoWeight() << std::endl;
                outfiles[precursor_isotope][fragment_isotope] << y_id.getContainer()[fragment_isotope].second
                                                              << "\t" << y_ion.getMonoWeight()
                                                              << "\t" << p.getMonoWeight() << std::endl;
            }
        }
    }
}


void sample_fragment_isotopic_distributions(std::string base_path, float max_mass, int num_samples, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium) {

    // create all output files and write header to each
    std::ofstream outfiles[max_depth][max_depth];
    for (int precursor_isotope = 1; precursor_isotope < max_depth; ++precursor_isotope) {
        for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
            std::string filename = "Precursor" + std::to_string(precursor_isotope) + "_" +
                                   "Fragment" + std::to_string(fragment_isotope) + ".tab";
            outfiles[precursor_isotope][fragment_isotope].open(base_path + filename);

            outfiles[precursor_isotope][fragment_isotope] << "probability" << "\tfrag.mass" << "\tprecursor.mass" << std::endl; //"\tfrag.a.mass" << "\tprecursor.a.mass" << std::endl;
        }
    }

    int max_length = max_mass/100;

    for (int peptide_length = 0; peptide_length <= max_length; ++peptide_length) {

        for (int sample = 0; sample < num_samples; ++sample) {

            OpenMS::AASequence random_sequence = create_random_peptide_sequence(peptide_length, num_sulfurs, num_c_sulfurs,
                                                                         num_selenium, num_c_selenium);

            if (random_sequence.size() > 0 && random_sequence.getMonoWeight() <= max_mass) {
                create_fragments(random_sequence, outfiles, num_sulfurs, num_c_sulfurs, num_selenium, num_c_selenium);
            }
        }

    }

    // close all output files
    for (int precursor_isotope = 1; precursor_isotope < max_depth; ++precursor_isotope) {
        for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
            outfiles[precursor_isotope][fragment_isotope].close();
        }
    }
}

int main(int argc, const char ** argv) {
    // Increase the maximum number of open files for this process. Was necessary for me.
    struct rlimit rlp;
    rlp.rlim_cur = 600;
    setrlimit(RLIMIT_NOFILE, &rlp);

    sample_fragment_isotopic_distributions(argv[1], atof(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));

    return 0;
}