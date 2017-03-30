//
// Created by Dennis Goldfarb on 1/3/17.
//

#include <iostream>
#include <fstream>
#include <random>

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

#include "FASTAParser.h"

static const ElementDB* elementDB = ElementDB::getInstance();



std::ofstream* openOutputFiles(std::string base_path, int max_depth)
{
    // create all output files and write header to each
    std::ofstream* outfiles = new std::ofstream[max_depth];
    for (int precursor_isotope = 0; precursor_isotope < max_depth; ++precursor_isotope)
    {
        std::string filename = "Precursor" + std::to_string(precursor_isotope) + ".tab";


        outfiles[precursor_isotope].open(base_path + filename);
        outfiles[precursor_isotope] << "probability" << "\tprecursor.mass" << std::endl;
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

void write_distribution(const OpenMS::AASequence &p, std::ofstream* outfiles, int max_depth, bool mono)
{
    OpenMS::EmpiricalFormula precursor_ef = p.getFormula();
    OpenMS::IsotopeDistribution precursor_id = precursor_ef.getIsotopeDistribution(0);

    for (int precursor_isotope = 0; precursor_isotope < max_depth; ++precursor_isotope)
    {
        double mass = mono ? precursor_ef.getMonoWeight() : precursor_ef.getAverageWeight();


        outfiles[precursor_isotope] << precursor_id.getContainer()[precursor_isotope].second << "\t" << mass << std::endl;
        //std::cout << precursor_id.getContainer()[precursor_isotope].second << "\t" << mass << std::endl;
    }
}

void sample_isotopic_distributions(std::string base_path, std::string fasta_path, float max_mass, int num_sulfurs, int max_depth, bool mono)
{
    std::ofstream* outfiles = openOutputFiles(base_path, max_depth);

    FASTAParser parser(fasta_path, max_mass, 1, 100);
    for (auto itr = parser.begin(); itr != parser.end(); ++itr)
    {
        if (num_sulfurs < 0 || itr->getFormula().getNumberOf(elementDB->getElement("Sulfur")) == num_sulfurs)
        {
            write_distribution(*itr, outfiles, max_depth, mono);
        }

    }

    closeOutputFiles(outfiles, max_depth);
}



void usage()
{
    std::cout << "GenerateTrainingData fasta_path out_path max_mass num_samples S max_depth mono" << std::endl;
    std::cout << "fasta_path: The path to the fasta file to train the splines on." << std::endl;
    std::cout << "out_path: The path to the directory that will store the training data, e.g. ~/data/" << std::endl;
    std::cout << "max_mass: maximum mass allowed for sampled peptides, e.g. 8500" << std::endl;
    std::cout << "S: number of sulfurs that should be in the fragment ion. Use -1 for all (e.g. 0,1,2..)" << std::endl;
    std::cout << "max_depth: The number of isotopes to generate training data for, e.g. 3 = M0,M1,M2" << std::endl;
    std::cout << "mono: should monoisotopic masses be used or average? 1=mono, 0=average" << std::endl;
    std::cout << std::endl;
}


int main(int argc, const char ** argv)
{
    if (argc != 7)
    {
        usage();
        return 0;
    }

    std::string fasta_path = argv[1];
    std::string out_path = argv[2];
    float max_mass = atof(argv[3]);
    int S = atoi(argv[4]);
    int max_depth = atoi(argv[5]);
    bool mono = strncmp(argv[6], "1", 1) == 0 ? true : false;

    std::cout << mono << std::endl;

    sample_isotopic_distributions(out_path, fasta_path, max_mass, S, max_depth, mono);

    return 0;
}