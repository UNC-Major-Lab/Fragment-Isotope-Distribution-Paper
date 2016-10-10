#include <iostream>
#include <random>
#include <numeric>

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include "Stats.h"

using namespace OpenMS;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

std::vector<double> sampleDecoy(int length)
{
    std::vector<double> probabilities;
    double sum = 0;
    for (int i = 0; i < length; ++i)
    {
        probabilities.push_back(dis(gen));
        sum+= probabilities[i];
    }
    for (int i = 0; i < length; ++i)
    {
        probabilities[i] /= sum;
    }

    return probabilities;
}

std::vector<double> sampleFromDistribution(std::vector<double>& probabilities)
{
    static int SAMPLE_SIZE = 1e3;

    std::vector<double> prefix_sum(probabilities.size());
    std::partial_sum(probabilities.begin(), probabilities.end(), prefix_sum.begin());

    std::vector<double> sample(probabilities.size(),0);

    for (int i = 0; i < SAMPLE_SIZE; ++i)
    {
        double ran = dis(gen);
        int index = std::upper_bound(prefix_sum.begin(), prefix_sum.end(), ran) - prefix_sum.begin();
        sample[index] += 1.0/SAMPLE_SIZE;
    }
    return sample;
}

std::vector<double> fillProbabilities(IsotopeDistribution& dist, UInt length)
{
    std::vector<double> probabilities;
    for (Size i = 0; i < dist.size(); ++i)
    {
        probabilities.push_back(dist.getContainer()[i].second);
    }

    for (Size i = dist.size(); i < length; ++i) {
        probabilities.push_back(0);
    }

    return probabilities;
}

std::vector<double> calculateScores(std::vector<double>& l, std::vector<double>& r)
{
    std::vector<double> result;
    result.push_back(Math::pearsonCorrelationCoefficient(l.begin(), l.end(), r.begin(), r.end()));
    result.push_back(Stats::totalVariationDistance(l.begin(), l.end(), r.begin(), r.end()));
    result.push_back(Stats::chiSquared(l.begin(), l.end(), r.begin(), r.end()));
    return result;
}

void testTheoreticalIon(AASequence& pep, AASequence& frag, EmpiricalFormula& precursor, EmpiricalFormula& fragment)
{
    static Size MAX_ISOTOPE = 3;
    static Size LENGTH = MAX_ISOTOPE+1;
    static Size num_processed = 0;

    int num_s = fragment.getNumberOf(ElementDB::getInstance()->getElement("Sulfur"));
    int num_cs = precursor.getNumberOf(ElementDB::getInstance()->getElement("Sulfur")) - num_s;

    std::vector<UInt> isolated_precursor_isotopes(1,0);
    for (UInt i = 1; i <= MAX_ISOTOPE; ++i) {
        isolated_precursor_isotopes.push_back(i);
        IsotopeDistribution frag_dist = fragment.getConditionalFragmentIsotopeDist(precursor, isolated_precursor_isotopes);

        IsotopeDistribution dist_averagine_precursor(i+1);
        dist_averagine_precursor.estimateFromPeptideWeight(fragment.getAverageWeight());
        dist_averagine_precursor.renormalize();
        IsotopeDistribution dist_exact_precursor = fragment.getIsotopeDistribution(i+1);

        std::vector<double> prob_averagine_precursor = fillProbabilities(dist_averagine_precursor, i+1);
        std::vector<double> prob_exact_precursor = fillProbabilities(dist_averagine_precursor, i+1);
        std::vector<double> prob_decoy = sampleDecoy(i+1);
        std::vector<double> prob_exact_fragment =  fillProbabilities(frag_dist, i+1);
        std::vector<double> prob_sampled_exact_fragment = sampleFromDistribution(prob_exact_fragment);

        std::vector<double> scores;

        scores = calculateScores(prob_averagine_precursor, prob_exact_fragment);

        //std::cout << scores[0] << "\t" << scores[1] << "\t" << scores[2] << "\t" << pep_mass << "\t"
        //          << frag_mass << "\t" << i << "\t" << num_s << "\t" << num_cs << "\t" << 0 << std::endl;

        scores = calculateScores(prob_exact_precursor, prob_exact_fragment);

        //std::cout << scores[0] << "\t" << scores[1] << "\t" << scores[2] << "\t" << pep_mass << "\t"
        //          << frag_mass << "\t" << i << "\t" << num_s << "\t" << num_cs << "\t" << 0 << std::endl;

        scores = calculateScores(prob_decoy, prob_exact_fragment);

        //std::cout << scores[0] << "\t" << scores[1] << "\t" << scores[2] << "\t" << pep_mass << "\t"
        //          << frag_mass << "\t" << i << "\t" << num_s << "\t" << num_cs << "\t" << 1 << std::endl;

        scores = calculateScores(prob_sampled_exact_fragment, prob_exact_fragment);

        //std::cout << scores[0] << "\t" << scores[1] << "\t" << scores[2] << "\t" << pep_mass << "\t"
        //          << frag_mass << "\t" << i << "\t" << num_s << "\t" << num_cs << "\t" << 0 << std::endl;
    }

    num_processed++;
    if (num_processed % 10000 == 0) {
        std::cout << num_processed << std::endl;
    }
}

void testTheoreticalPeptide(AASequence& pep)
{
    EmpiricalFormula precursor = pep.getFormula();
    EmpiricalFormula fragment;
    for (Size i = 1; i < pep.size(); i++)
    {
        AASequence frag = pep.getPrefix(i);

        fragment = frag.getFormula(Residue::ResidueType::BIon);
        testTheoreticalIon(pep, frag, precursor, fragment);

        fragment = pep.getPrefix(i).getFormula(Residue::ResidueType::YIon);
        testTheoreticalIon(pep, frag, precursor, fragment);
    }
}

void testTheoreticalProtein(FASTAFile::FASTAEntry& protein, EnzymaticDigestion& digestor)
{
    static Size MIN_PEPTIDE_LENGTH = 5;
    static Size MAX_PEPTIDE_LENGTH = 80;

    std::vector<AASequence> peptides;
    digestor.digest(AASequence::fromString(protein.sequence), peptides);
    for (Size j = 0; j < peptides.size(); ++j)
    {
        if (peptides[j].size() >= MIN_PEPTIDE_LENGTH && peptides[j].size() <= MAX_PEPTIDE_LENGTH)
        {
            testTheoreticalPeptide(peptides[j]);
        }
    }
}

void testTheoreticalPeptides(std::string fasta_path)
{
    std::vector<FASTAFile::FASTAEntry> proteins;
    FASTAFile().load(fasta_path, proteins);

    EnzymaticDigestion digestor; // default parameters are fully tryptic with 0 missed cleavages

    for (Size i = 0; i < proteins.size(); ++i)
    {
        testTheoreticalProtein(proteins[i], digestor);
    }
}

int main(int argc, char * argv[])
{
    testTheoreticalPeptides(argv[1]);
    return 0;
}