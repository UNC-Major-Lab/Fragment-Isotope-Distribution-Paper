#include <iostream>
#include <fstream>
#include <random>
#include <numeric>
#include <string>
#include <set>
#include <functional>
#include <string>

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include "Stats.h"

using namespace OpenMS;

static const ElementDB* elementDB = ElementDB::getInstance();
static const IsotopeSplineDB* isotopeDB = IsotopeSplineDB::getInstance();

static Size MAX_ISOTOPE = 4;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

std::set<AASequence> uniquePeptides;

std::map<std::string, std::pair<std::vector<double>, std::vector<double> > > precursor_method2val;
std::map<std::string, std::map<std::string, std::pair<std::vector<double>, std::vector<double> > > > fragment_method2iso2val;


bool isValidPeptide(AASequence& pep) {
    String p = pep.toString();
    if (p.hasSubstring("U") || p.hasSubstring("B") || p.hasSubstring("Z") || p.hasSubstring("J") || p.hasSubstring("X"))
    {
        return false;
    }
    return true;
}

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

std::vector<double> calculateResiduals(std::vector<double>& l, std::vector<double>& r)
{
    std::vector<double> result;
    for (int i = 0; i < l.size(); ++i)
    {
        result.push_back(r[i]-l[i]);
    }
    return result;
}

void testTheoreticalIsolation(EmpiricalFormula& precursor, EmpiricalFormula& fragment, std::set<UInt>& isolated_precursor_isotopes,
                              double pep_mass, double frag_mass, int num_s_prec, int num_s_frag, UInt depth, std::string label)
{
    IsotopeDistribution exact_fragment_dist = fragment.getConditionalFragmentIsotopeDist(precursor, isolated_precursor_isotopes);

    IsotopeDistribution approx_precursor_dist(depth);
    approx_precursor_dist.estimateFromPeptideWeight(frag_mass);
    approx_precursor_dist.renormalize();

    IsotopeDistribution approx_fragment_dist(depth);
    approx_fragment_dist.estimateForFragmentFromPeptideWeight(pep_mass, frag_mass, isolated_precursor_isotopes);
    approx_fragment_dist.renormalize();

    IsotopeDistribution approx_fragment_S_dist(depth);
    approx_fragment_S_dist.estimateForFragmentFromPeptideWeightAndS(pep_mass, num_s_prec, frag_mass, num_s_frag, isolated_precursor_isotopes);
    approx_fragment_S_dist.renormalize();

    IsotopeDistribution approx_fragment_spline_dist = isotopeDB->estimateForFragmentFromPeptideWeight(pep_mass, frag_mass, isolated_precursor_isotopes);
    approx_fragment_spline_dist.renormalize();

    IsotopeDistribution approx_fragment_splineS_dist = isotopeDB->estimateForFragmentFromPeptideWeightAndS(pep_mass, num_s_prec, frag_mass, num_s_frag, isolated_precursor_isotopes);
    approx_fragment_splineS_dist.renormalize();


    std::vector<double> exact_fragment_prob =  fillProbabilities(exact_fragment_dist, depth);
    std::vector<double> approx_precursor_prob = fillProbabilities(approx_precursor_dist, depth);
    std::vector<double> approx_fragment_prob = fillProbabilities(approx_fragment_dist, depth);
    std::vector<double> approx_fragment_S_prob = fillProbabilities(approx_fragment_S_dist, depth);
    std::vector<double> approx_fragment_spline_prob = fillProbabilities(approx_fragment_spline_dist, depth);
    std::vector<double> approx_fragment_splineS_prob = fillProbabilities(approx_fragment_spline_dist, depth);

    //std::vector<double> decoy_prob = sampleDecoy(i+1);
    //std::vector<double> sampled_exact_fragment_prob = sampleFromDistribution(exact_fragment_prob);

    std::vector<double> scores;

    //scores = calculateScores(exact_fragment_prob, approx_precursor_prob);
    //out_scores << scores[2] << "\t" << label << "\t" << "p" << std::endl;

    scores = calculateScores(exact_fragment_prob, approx_fragment_prob);
    fragment_method2iso2val["Averagine"][label].first.push_back(scores[2]);

    scores = calculateScores(exact_fragment_prob, approx_fragment_S_prob);
    fragment_method2iso2val["Sulfur-specific Averagine"][label].first.push_back(scores[2]);

    scores = calculateScores(exact_fragment_prob, approx_fragment_spline_prob);
    fragment_method2iso2val["Spline"][label].first.push_back(scores[2]);

    scores = calculateScores(exact_fragment_prob, approx_fragment_splineS_prob);
    fragment_method2iso2val["Sulfur-specific spline"][label].first.push_back(scores[2]);

    //Residuals
    //scores = calculateResiduals(exact_fragment_prob, approx_precursor_prob);
    //for (int i = 0; i < scores.size(); ++i) out_residual << scores[i] << "\t" << label << "\t" << "p" << std::endl;

    scores = calculateResiduals(exact_fragment_prob, approx_fragment_prob);
    for (int i = 0; i < scores.size(); ++i) fragment_method2iso2val["Averagine"][label].second.push_back(scores[i]);

    scores = calculateResiduals(exact_fragment_prob, approx_fragment_S_prob);
    for (int i = 0; i < scores.size(); ++i) fragment_method2iso2val["Sulfur-specific Averagine"][label].second.push_back(scores[i]);

    scores = calculateResiduals(exact_fragment_prob, approx_fragment_spline_prob);
    for (int i = 0; i < scores.size(); ++i) fragment_method2iso2val["Spline"][label].second.push_back(scores[i]);

    scores = calculateResiduals(exact_fragment_prob, approx_fragment_splineS_prob);
    for (int i = 0; i < scores.size(); ++i) fragment_method2iso2val["Sulfur-specific spline"][label].second.push_back(scores[i]);
}

void testTheoreticalIon(AASequence& pep, AASequence& frag, EmpiricalFormula& precursor, EmpiricalFormula& fragment)
{


    int num_s_frag = fragment.getNumberOf(ElementDB::getInstance()->getElement("Sulfur"));
    int num_s_prec = precursor.getNumberOf(ElementDB::getInstance()->getElement("Sulfur"));

    double pep_mass = precursor.getAverageWeight();
    double frag_mass = fragment.getAverageWeight();

    std::set<UInt> isolated_precursor_isotopes;
    for (UInt i = 1; i <= MAX_ISOTOPE; ++i) {
        isolated_precursor_isotopes.insert(i);
        std::string label = "0-"+std::to_string(i);
        testTheoreticalIsolation(precursor, fragment, isolated_precursor_isotopes, pep_mass, frag_mass, num_s_prec, num_s_frag, i+1, label);
    }

    for (UInt i = 1; i <= MAX_ISOTOPE; ++i) {
        isolated_precursor_isotopes.clear();
        isolated_precursor_isotopes.insert(i);
        testTheoreticalIsolation(precursor, fragment, isolated_precursor_isotopes, pep_mass, frag_mass, num_s_prec, num_s_frag, i+1, std::to_string(i));
    }
}

void testTheoreticalPeptideDistribution(EmpiricalFormula &p)
{
    UInt depth = 11;
    IsotopeDistribution exact, averagine(depth), spline(depth), averagineS(depth), splineS(depth);

    int num_S = p.getNumberOf(elementDB->getElement("Sulfur"));

    double average_weight = p.getAverageWeight();
    exact = p.getIsotopeDistribution(depth);
    averagine.estimateFromPeptideWeight(average_weight);
    averagineS.estimateFromPeptideWeightAndS(average_weight, num_S);
    //averagine.estimateFromWeightAndComp(average_weight, 4.86151, 7.68282, 1.3005, 1.56299, 0.047074, 0);
    spline = isotopeDB->estimateFromPeptideWeight(average_weight, depth);
    splineS = isotopeDB->estimateFromPeptideWeightAndS(average_weight, num_S, depth);

    std::vector<double> exact_prob =  fillProbabilities(exact, depth);
    std::vector<double> averagine_prob = fillProbabilities(averagine, depth);
    std::vector<double> averagineS_prob = fillProbabilities(averagineS, depth);
    std::vector<double> spline_prob =  fillProbabilities(spline, depth);
    std::vector<double> splineS_prob =  fillProbabilities(splineS, depth);

    std::vector<double> scores;
    scores = calculateScores(exact_prob, averagine_prob);
    precursor_method2val["Averagine"].first.push_back(scores[2]);
    scores = calculateScores(exact_prob, averagineS_prob);
    precursor_method2val["Sulfur-specific averagine"].first.push_back(scores[2]);
    scores = calculateScores(exact_prob, spline_prob);
    precursor_method2val["Spline"].first.push_back(scores[2]);
    scores = calculateScores(exact_prob, splineS_prob);
    precursor_method2val["Sulfur-specific spline"].first.push_back(scores[2]);
    /*scores = calculateScores(averagine_prob, spline_prob);
    out_scores << scores[2] << "\t" << average_weight << "\t" << "averagine vs spline" << std::endl;
    */


    scores = calculateResiduals(exact_prob, averagine_prob);
    for (int i = 0; i < scores.size(); ++i) precursor_method2val["Averagine"].second.push_back(scores[i]);
    scores = calculateResiduals(exact_prob, averagineS_prob);
    for (int i = 0; i < scores.size(); ++i) precursor_method2val["Sulfur-specific averagine"].second.push_back(scores[i]);
    scores = calculateResiduals(exact_prob, spline_prob);
    for (int i = 0; i < scores.size(); ++i) precursor_method2val["Spline"].second.push_back(scores[i]);
    scores = calculateResiduals(exact_prob, splineS_prob);
    for (int i = 0; i < scores.size(); ++i) precursor_method2val["Sulfur-specific spline"].second.push_back(scores[i]);
    /*scores = calculateResiduals(averagine_prob, spline_prob);
    for (int i = 0; i < scores.size(); ++i) out_residual << scores[i] << "\t" << "averagine vs spline" << std::endl;
    */
}

void testTheoreticalPeptide(AASequence& pep, bool doFragments)
{
    EmpiricalFormula precursor = pep.getFormula();
    EmpiricalFormula fragment;

    if (doFragments)
    {
        for (Size i = 1; i < pep.size(); i++)
        {
            AASequence frag = pep.getPrefix(i);

            fragment = frag.getFormula(Residue::ResidueType::BIon);
            testTheoreticalIon(pep, frag, precursor, fragment);

            fragment = pep.getPrefix(i).getFormula(Residue::ResidueType::YIon);
            testTheoreticalIon(pep, frag, precursor, fragment);
        }
    }
    else
    {
        testTheoreticalPeptideDistribution(precursor);
    }
}

void testTheoreticalProtein(FASTAFile::FASTAEntry& protein, EnzymaticDigestion& digestor, bool doFragments)
{
    static Size MIN_PEPTIDE_LENGTH = 5;
    static Size MAX_PEPTIDE_LENGTH = 80;

    std::vector<AASequence> peptides;
    digestor.digest(AASequence::fromString(protein.sequence), peptides);
    for (Size j = 0; j < peptides.size(); ++j)
    {
        if (peptides[j].size() >= MIN_PEPTIDE_LENGTH && peptides[j].size() <= MAX_PEPTIDE_LENGTH
            && isValidPeptide(peptides[j]) && uniquePeptides.find(peptides[j]) == uniquePeptides.end())
        {
            uniquePeptides.insert(peptides[j]);
    		testTheoreticalPeptide(peptides[j], doFragments);
        }
    }
}

void testTheoreticalPeptides(std::string fasta_path, int job_id, int num_jobs, bool doFragments)
{
    std::vector<FASTAFile::FASTAEntry> proteins;
    FASTAFile().load(fasta_path, proteins);

    EnzymaticDigestion digestor; // default parameters are fully tryptic with 0 missed cleavages

    for (Size i = job_id; i < proteins.size(); i+=num_jobs)
    {
        testTheoreticalProtein(proteins[i], digestor, doFragments);
    }
}

void writeResults(std::string path_residual, std::string path_chisquared, std::string path_stats, bool doFragments)
{
    std::ofstream out_residual(path_residual);
    std::ofstream out_scores(path_chisquared);
    std::ofstream out_stats(path_stats);

    if (doFragments)
    {
        std::map<std::string, std::map<std::string, std::map<int, int> > > fragment_method2iso2bin2count;
        for (auto const &method_itr : fragment_method2iso2val)
        {
            std::string const &key = method_itr.first;
            for (auto const &iso_itr : fragment_method2iso2val[key])
            {
                std::string const &iso = iso_itr.first;
                std::vector<double> chi = iso_itr.second.first;
                std::vector<double> res = iso_itr.second.second;
                std::sort(chi.begin(), chi.end());
                std::sort(res.begin(), res.end());

                double mean_chi = std::accumulate(chi.begin(), chi.end(), 0.0, std::plus<double>()) / chi.size();
                for (int i = 0; i < res.size(); ++i) res[i] = std::abs(res[i]);
                double mean_res = std::accumulate(res.begin(), res.end(), 0.0, std::plus<double>()) / res.size();

                double median_chi = chi[chi.size()/2];
                double median_res = res[res.size()/2];

                double q1_chi = chi[chi.size()/4];
                double q1_res = res[res.size()/4];

                double q3_chi = chi[3*chi.size()/4];
                double q3_res = res[3*res.size()/4];

                double min_chi = chi[0];
                double min_res = res[0];

                double max_chi = chi[chi.size()-1];
                double max_res = res[res.size()-1];

                out_stats << mean_chi << "\t" << min_chi << "\t" << q1_chi << "\t" << median_chi << "\t" << q3_chi << "\t" << max_chi << "\t" << iso << "\t" << key << std::endl;
                out_stats << mean_res << "\t" << min_res << "\t" << q1_res << "\t" << median_res << "\t" << q3_res << "\t" << max_res << "\t" << iso << "\t" << key << std::endl;

            }
        }
    } else
    {
        std::map<std::string, std::map<int, int> > precursor_method2bin2count;

        for (auto const &method_itr : precursor_method2val)
        {
            std::string const &key = method_itr.first;
            std::vector<double> chi = method_itr.second.first;
            std::vector<double> res = method_itr.second.second;
            double mean_chi = std::accumulate(chi.begin(), chi.end(), 0.0, std::plus<double>()) / chi.size();

            for (int i = 0; i < res.size(); ++i) res[i] = std::abs(res[i]);
            double mean_res = std::accumulate(res.begin(), res.end(), 0.0, std::plus<double>()) / res.size();

            std::sort(chi.begin(), chi.end());
            std::sort(res.begin(), res.end());

            double median_chi = chi[chi.size()/2];
            double median_res = res[res.size()/2];

            double q1_chi = chi[chi.size()/4];
            double q1_res = res[res.size()/4];

            double q3_chi = chi[3*chi.size()/4];
            double q3_res = res[3*res.size()/4];

            double min_chi = chi[0];
            double min_res = res[0];

            double max_chi = chi[chi.size()-1];
            double max_res = res[res.size()-1];

            out_stats << mean_chi << "\t" << min_chi << "\t" << q1_chi << "\t" << median_chi << "\t" << q3_chi << "\t" << max_chi << "\t" << key << std::endl;
            out_stats << mean_res << "\t" << min_res << "\t" << q1_res << "\t" << median_res << "\t" << q3_res << "\t" << max_res << "\t" << key << std::endl;

        }
    }

    out_residual.close();
    out_scores.close();
    out_stats.close();
}

void init(bool doFragments)
{

    if (doFragments)
    {
        for (UInt i = 1; i <= MAX_ISOTOPE; ++i) {
            std::string label = "0-" + std::to_string(i);
            fragment_method2iso2val["Averagine"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
            fragment_method2iso2val["Sulfur-specific averagine"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
            fragment_method2iso2val["Spline"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
            fragment_method2iso2val["Sulfurs-specific spline"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
            label = std::to_string(i);
            fragment_method2iso2val["Averagine"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
            fragment_method2iso2val["Sulfur-specific averagine"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
            fragment_method2iso2val["Spline"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
            fragment_method2iso2val["Sulfurs-specific spline"][label] = std::make_pair(std::vector<double>(), std::vector<double>());
        }
    } else
    {
        precursor_method2val["Averagine"] = std::make_pair(std::vector<double>(), std::vector<double>());
        precursor_method2val["Sulfur-specific averagine"] = std::make_pair(std::vector<double>(), std::vector<double>());
        precursor_method2val["Spline"] = std::make_pair(std::vector<double>(), std::vector<double>());
        precursor_method2val["Sulfurs-specific spline"] = std::make_pair(std::vector<double>(), std::vector<double>());
    }


}

void usage()
{
    std::cout << "CompareToTheoretical fasta_path job_id num_jobs do_frag residual_file score_file stats_file" << std::endl;
}

int main(int argc, char * argv[])
{
    if (argc != 7)
    {
        usage();
    }

    init(atoi(argv[4]));

    testTheoreticalPeptides(argv[1], atoi(argv[2])-1, atoi(argv[3]), atoi(argv[4]));

    writeResults(argv[5], argv[6], argv[7], atoi(argv[4]));



    return 0;
}
