//
// Created by Dennis Goldfarb on 1/3/17.
//

#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

#include "FASTAParser.h"

static const ElementDB* elementDB = ElementDB::getInstance();
static const OpenMS::ResidueDB* residueDB = OpenMS::ResidueDB::getInstance();

static std::string AMINO_ACIDS = "ADEFGHIKLNPQRSTVWYCM";
static std::string AMINO_ACIDS_NO_SULFUR = "ADEFGHIKLNPQRSTVWY";
static std::string AMINO_ACIDS_SULFUR = "CM";

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_AA(0,1);
std::uniform_int_distribution<> dis_SULFUR(0,1);

std::ofstream* openOutputFiles(std::string base_path, int max_depth, bool write_sulfur)
{
    // create all output files and write header to each
    std::ofstream* outfiles = new std::ofstream[max_depth];
    for (int precursor_isotope = 0; precursor_isotope < max_depth; ++precursor_isotope)
    {
        std::string filename = "Precursor" + std::to_string(precursor_isotope) + ".tab";


        outfiles[precursor_isotope].open(base_path + filename);
        if (write_sulfur) {
            outfiles[precursor_isotope] << "probability" << "\tprecursor.mass" << "\tsulfur" << std::endl;
        } else {
            outfiles[precursor_isotope] << "probability" << "\tprecursor.mass" << std::endl;
        }
    }

    return outfiles;
}

void closeOutputFiles(std::ofstream* outfiles, int max_depth)
{
    // close all output files
    for (int precursor_isotope = 0; precursor_isotope < max_depth; ++precursor_isotope)
    {
        outfiles[precursor_isotope].close();
    }
    delete[] outfiles;
}

void write_distribution(const OpenMS::AASequence &p, std::ofstream* outfiles, int max_depth, bool mono, bool write_sulfur, bool exact)
{
    OpenMS::EmpiricalFormula precursor_ef = p.getFormula();
    double mass = mono ? precursor_ef.getMonoWeight() : precursor_ef.getAverageWeight();

    OpenMS::IsotopeDistribution precursor_id;
    if (exact) {
        precursor_id = precursor_ef.getIsotopeDistribution(0);
    } else {
        precursor_id.estimateFromPeptideWeight(mass);
    }

    for (int precursor_isotope = 0; precursor_isotope < max_depth && precursor_isotope < precursor_id.size(); ++precursor_isotope)
    {
        if (write_sulfur) {
            int num_sulfur = precursor_ef.getNumberOf(elementDB->getElement("Sulfur"));
            outfiles[precursor_isotope] << precursor_id.getContainer()[precursor_isotope].second << "\t" << mass << "\t" << num_sulfur << std::endl;
        } else {
            outfiles[precursor_isotope] << precursor_id.getContainer()[precursor_isotope].second << "\t" << mass << std::endl;
        }
    }
}

void proteome_isotopic_distributions(std::string base_path, std::string fasta_path, float max_mass, int max_depth, bool mono)
{
    std::ofstream* outfiles = openOutputFiles(base_path, max_depth, true);
    std::ofstream* outfiles_averagine = openOutputFiles(base_path+"averagine/", max_depth, true);

    FASTAParser parser(fasta_path, max_mass, 1, 150);
    for (auto itr = parser.begin(); itr != parser.end(); ++itr)
    {
        write_distribution(*itr, outfiles, max_depth, mono, true, true);
        write_distribution(*itr, outfiles, max_depth, mono, true, false);
    }

    closeOutputFiles(outfiles, max_depth);
    closeOutputFiles(outfiles_averagine, max_depth);
}

OpenMS::AASequence create_random_peptide_sequence(int peptide_length, std::vector<double> aa2prob, int num_sulfurs)
{
    OpenMS::AASequence random_peptide;

    if (num_sulfurs < 0)
    {
        for (int aa_index = 0; aa_index < peptide_length; ++aa_index)
        {
            double rand = dis_AA(gen);
            int index = std::lower_bound(aa2prob.begin(), aa2prob.end(), rand) - aa2prob.begin() - 1;
            random_peptide += residueDB->getResidue(AMINO_ACIDS[index]);
        }
    }
    else
    {
        // for insertion of sulfur containing amino acids
        for (int i = 0; i < num_sulfurs; ++i)
        {
            random_peptide += residueDB->getResidue(AMINO_ACIDS_SULFUR[dis_SULFUR(gen)]);
        }

        // random amino acid insertion (non Sulfur and Selenium amino acids)
        for (int aa_index = 0; aa_index < peptide_length; ++aa_index)
        {
            double rand = dis_AA(gen);
            int index = std::lower_bound(aa2prob.begin(), aa2prob.end(), rand) - aa2prob.begin() - 1;
            random_peptide += residueDB->getResidue(AMINO_ACIDS_NO_SULFUR[index]);
        }
    }

    return random_peptide;
}


std::map<char, double> getAAProbabilities(std::string fasta_path, bool sulfur)
{
    std::map<char, double> aa2prob;
    for (char aa : AMINO_ACIDS)
    {
        aa2prob[aa] = 0.0;
    }

    std::vector<FASTAFile::FASTAEntry> proteins;
    FASTAFile().load(fasta_path, proteins);

    int count = 0;
    for (Size i = 0; i < proteins.size(); ++i)
    {
        for (int j = 0; j < proteins[i].sequence.size(); ++j)
        {
            char aa = proteins[i].sequence[j];
            if ((sulfur && AMINO_ACIDS.find(aa) != -1) || (!sulfur && AMINO_ACIDS_NO_SULFUR.find(aa) != -1))
            {
                aa2prob[aa]++;
                count++;
            }
        }
    }

    for (auto &aa : aa2prob)
    {
        aa.second /= count;
    }

    return aa2prob;
}

std::vector<double> calcPrefixSum(std::map<char, double> aa2prob, bool sulfur)
{
    std::string AAs = sulfur ? AMINO_ACIDS : AMINO_ACIDS_NO_SULFUR;
    std::vector<double> prefixSum;

    prefixSum.push_back(0);

    for (int i = 0; i < AAs.size(); ++i)
    {
        prefixSum.push_back(aa2prob[AAs[i]] + prefixSum[i]);
    }

    return prefixSum;
}

void sample_isotopic_distributions(std::string base_path, std::string fasta_path, float max_mass, int num_sulfurs, int num_samples, int max_depth, bool mono)
{

    std::vector<double> aa2prob = calcPrefixSum(getAAProbabilities(fasta_path, num_sulfurs == -1), num_sulfurs == -1);

    std::ofstream* outfiles = openOutputFiles(base_path, max_depth, false);

    int max_length = max_mass/100;


    for (int peptide_length = 1; peptide_length <= max_length; ++peptide_length)
    {
        for (int sample = 0; sample < num_samples; ++sample)
        {
            OpenMS::AASequence random_sequence = create_random_peptide_sequence(peptide_length, aa2prob, num_sulfurs);

            if (random_sequence.size() > 0 && random_sequence.getMonoWeight() <= max_mass)
            {
                write_distribution(random_sequence, outfiles, max_depth, mono, false, true);
            }
        }
    }

    closeOutputFiles(outfiles, max_depth);
}






void usage()
{
    std::cout << "GenerateTrainingData fasta_path out_path max_mass S num_samples max_depth mono" << std::endl;
    std::cout << "fasta_path: The path to the fasta file to train the splines on." << std::endl;
    std::cout << "out_path: The path to the directory that will store the training data, e.g. ~/data/" << std::endl;
    std::cout << "max_mass: maximum mass allowed for sampled peptides, e.g. 8500" << std::endl;
    std::cout << "max_depth: The number of isotopes to generate training data for, e.g. 3 = M0,M1,M2" << std::endl;
    std::cout << "mono: should monoisotopic masses be used or average? 1=mono, 0=average" << std::endl;
    std::cout << "S: number of sulfurs that should be in the fragment ion. Use -1 for all (e.g. 0,1,2..)" << std::endl;
    std::cout << "num_samples: number of random peptides to make for each peptide length" << std::endl;

    std::cout << std::endl;
}


int main(int argc, const char ** argv)
{
    if (argc != 8 && argc != 6)
    {
        usage();
        return 0;
    }

    std::string fasta_path = argv[1];
    std::string out_path = argv[2];
    float max_mass = atof(argv[3]);
    int max_depth = atoi(argv[4]);
    bool mono = strncmp(argv[5], "1", 1) == 0 ? true : false;

    if (argc == 8) {
        int S = atoi(argv[6]);
        int num_samples = atoi(argv[7]);
        sample_isotopic_distributions(out_path, fasta_path, max_mass, S, num_samples, max_depth, mono);
    } else {
        proteome_isotopic_distributions(out_path, fasta_path, max_mass, max_depth, mono);
    }

    return 0;
}