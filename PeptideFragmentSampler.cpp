//
// Created by Dennis Goldfarb on 9/28/16.
//

#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <algorithm>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

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

int max_depth;

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

OpenMS::IsotopeDistribution getIsotopeDistribution(OpenMS::EmpiricalFormula &formula, OpenMS::UInt max_depth) const
{
    OpenMS::IsotopeDistribution result(max_depth);
    for (auto it = formula.begin(); it != formula.end(); ++it)
    {
        OpenMS::IsotopeDistribution tmp = it->first->getIsotopeDistribution();
        tmp.setMaxIsotope(max_depth);
        result += tmp * it->second;
    }
    return result;
}

OpenMS::IsotopeDistribution getFragmentDistribution(OpenMS::EmpiricalFormula &precursor_ef, OpenMS::EmpiricalFormula &fragment_ef, std::vector<OpenMS::UInt> &isolated_isotopes)
{
    // A fragment's isotopes can only be as high as the largest isolated precursor isotope.
    OpenMS::UInt max_depth = *std::max_element(isolated_isotopes.begin(), isolated_isotopes.end())+1;

    // Treat *this as the fragment molecule
    OpenMS::EmpiricalFormula complementary_fragment = precursor_ef-fragment_ef;

    OpenMS::IsotopeDistribution fragment_isotope_dist = getIsotopeDistribution(fragment_ef, max_depth);
    OpenMS::IsotopeDistribution comp_fragment_isotope_dist = getIsotopeDistribution(complementary_fragment, max_depth);

    OpenMS::IsotopeDistribution result;
    result.calcFragmentIsotopeDist(fragment_isotope_dist, comp_fragment_isotope_dist, isolated_isotopes);

    return result;
}

void create_fragments(OpenMS::AASequence &p, std::ofstream** outfiles, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium) {
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

            //OpenMS::IsotopeDistribution b_id = b_ion.getConditionalFragmentIsotopeDist(precursor_ef, isolated_isotopes);
            //OpenMS::IsotopeDistribution y_id = y_ion.getConditionalFragmentIsotopeDist(precursor_ef, isolated_isotopes);
            OpenMS::IsotopeDistribution b_id = getFragmentDistribution(precursor_ef, b_ion, isolated_isotopes);
            OpenMS::IsotopeDistribution y_id = getFragmentDistribution(precursor_ef, y_ion, isolated_isotopes);

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


void sample_fragment_isotopic_distributions(std::string base_path, float max_mass, int num_samples, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium, bool append) {

    // create all output files and write header to each
    std::ofstream** outfiles = new std::ofstream*[max_depth];
    for (int precursor_isotope = 0; precursor_isotope < max_depth; ++precursor_isotope) {
        outfiles[precursor_isotope] = new std::ofstream[max_depth];
        for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
            std::string filename = "Precursor" + std::to_string(precursor_isotope) + "_" +
                                   "Fragment" + std::to_string(fragment_isotope) + ".tab";

            if (append)
            {
                outfiles[precursor_isotope][fragment_isotope].open(base_path + filename, std::ofstream::out | std::ofstream::app);
            } else
            {
                outfiles[precursor_isotope][fragment_isotope].open(base_path + filename);
            }

            // Only add the header if we're creating a new file
            // This happens if we're not appending, or if we're appending and it's our first time in the loop
            if (!append) {
                outfiles[precursor_isotope][fragment_isotope] << "probability" << "\tfrag.mass" << "\tprecursor.mass"
                                                              << std::endl; //"\tfrag.a.mass" << "\tprecursor.a.mass" << std::endl;
            }
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
    for (int precursor_isotope = 0; precursor_isotope < max_depth; ++precursor_isotope) {
        for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
            outfiles[precursor_isotope][fragment_isotope].close();
        }
        delete[] outfiles[precursor_isotope];
    }
    delete[] outfiles;
}

void sample_average_fragment_isotopic_distribution(std::string distribution_path, std::string base_path, float max_mass, float min_percentage)
{
    std::map<std::pair<int,int>, int> sulfurs2count;

    std::ifstream sulfur_dist_in(distribution_path);
    std::string input;
    int S, CS, count, i=0;
    while (sulfur_dist_in >> input)
    {
        if (i == 0) S = atoi(input.c_str());
        else if (i == 1) CS = atoi(input.c_str());
        else {
            count = atoi(input.c_str());
            sulfurs2count[std::make_pair(S,CS)] = count;
        }
        i = (i+1) % 3;
    }

    int max_count = 0;
    for (auto itr : sulfurs2count)
    {
        max_count = std::max(max_count, itr.second);
    }

    bool append = false;
    for (auto itr : sulfurs2count)
    {
        double percentage = (double) itr.second / max_count;
        if (percentage >= min_percentage) {
            int num_samples = std::floor(percentage / min_percentage);
            sample_fragment_isotopic_distributions(base_path, max_mass, num_samples, itr.first.first, itr.first.second, 0, 0, append);
            append = false;
        }
    }
}

void usage()
{
    std::cout << "PeptideFragmentSampler out_path max_mass num_samples S CS Se CSe max_isotope" << std::endl;
    std::cout << "out_path: The path to the directory that will store the training data, e.g. ~/data/" << std::endl;
    std::cout << "max_mass: maximum mass allowed for sampled peptides, e.g. 8500" << std::endl;
    std::cout << "num_samples: number of random peptides to generate for each peptide length, e.g 100" << std::endl;
    std::cout << "S: number of sulfurs that should be in the fragment ion, e.g. 0" << std::endl;
    std::cout << "CS: number of sulfurs that should be in the complementary fragment ion, e.g. 0" << std::endl;
    std::cout << "Se: number of seleniums that should be in the fragment ion, e.g. 0" << std::endl;
    std::cout << "CSe: number of seleniums that should be in the complementary fragment ion, e.g. 0" << std::endl;
    std::cout << "max_isotope: The maximum isotope to generate training data for, e.g. 5" << std::endl;
    std::cout << std::endl;

    std::cout << "PeptideFragmentSampler sulfur_dist_path out_path max_mass min_percentage max_isotope" << std::endl;
    std::cout << "sulfur_dist_path: file path to the results of GetSulfurDistribution, e.g. ~/data/sulfur_distribution.tab" << std::endl;
    std::cout << "out_path: The path to the directory that will store the training data, e.g. ~/data/" << std::endl;
    std::cout << "max_mass: maximum mass allowed for sampled peptides, e.g. 8500" << std::endl;
    std::cout << "min_percentage: the min abundance of a sulfur distribution necessary to be included in the training data (relative to most abundant case), e.g. .001" << std::endl;
    std::cout << "max_isotope: The maximum isotope to generate training data for, e.g. 5" << std::endl;
}

int main(int argc, const char ** argv) {

    if (argc != 9 && argc != 6)
    {
        usage();
    }

    // Increase the maximum number of open files for this process. Was necessary for me.
    struct rlimit rlp;
    rlp.rlim_cur = 600;
    setrlimit(RLIMIT_NOFILE, &rlp);

    if (argc == 9) {
        std::string out_path = argv[1];
        float max_mass = atof(argv[2]);
        int num_samples = atoi(argv[3]);
        int S = atoi(argv[4]);
        int CS = atoi(argv[5]);
        int Se = atoi(argv[6]);
        int CSe = atoi(argv[7]);

        max_depth = atoi(argv[8]) + 1;

        sample_fragment_isotopic_distributions(out_path, max_mass, num_samples, S, CS, Se, CSe, false);
    }
    else if (argc == 6)
    {
        std::string dist_path = argv[1];
        std::string out_path = argv[2];
        float max_mass = atof(argv[3]);
        float min_percentage = atof(argv[4]);

        max_depth = atoi(argv[5]) + 1;

        sample_average_fragment_isotopic_distribution(dist_path, out_path, max_mass, min_percentage);
    }

    return 0;
}