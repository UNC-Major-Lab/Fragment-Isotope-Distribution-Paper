//
// Created by Dennis Goldfarb on 11/14/16.
//
#include <iostream>
#include <fstream>
#include <string>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>

#include "Ion.h"
#include "SpectrumUtilities.h"
#include "Stats.h"

static const OpenMS::IsotopeSplineDB* isotopeDB = OpenMS::IsotopeSplineDB::getInstance();

void normalizeDist(std::vector<std::pair<double, double> > &dist)
{
    double maxProb = 0;
    for (auto itr = dist.begin(); itr != dist.end(); ++itr)
    {
        maxProb = std::max(maxProb, itr->second);
    }
    for (int i = 0; i < dist.size(); ++i)
    {
        dist[i].second /= maxProb;
    }
}

void outputDist(std::ofstream &out, std::vector<std::pair<double, double> > &dist, std::string ion_name, int ionIndex,
                std::string isotope_range, std::string name)
{


    for (int i = 0; i < dist.size(); ++i)
    {
        out << isotope_range << "\t" << ionIndex << "\t" << ion_name << "\t" << dist[i].first << "\t"
                 << dist[i].second << "\t" << name << std::endl;
    }
}

void outputScores(std::ofstream &out, std::string ion_name, int ionIndex, std::string isotope_range, std::string name, double score, double x, double y)
{
    out << isotope_range << "\t" << ionIndex << "\t" << ion_name << "\t" << name << "\t" << score << "\t" << x << "\t" << y << std::endl;
}


void usage()
{
    std::cout << "Usage: " << std::endl;
}

int main(int argc, char * argv[])
{

    OpenMS::MzMLFile mzMLDataFileProfile, mzMLDataFileCentroid;
    OpenMS::MSExperiment<OpenMS::Peak1D> msExperimentProfile, msExperimentCentroid;
    try {
        mzMLDataFileProfile.load(argv[1], msExperimentProfile);
        mzMLDataFileCentroid.load(argv[2], msExperimentCentroid);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 1;
    }

    std::ofstream out(argv[3]);
    std::ofstream calc_out(argv[4]);
    std::ofstream scores_out(argv[5]);

    out << "isotope.range" << "\t" << "ion.index" << "\t" << "ion.name" << "\t" << "mz" << "\t" << "int" << std::endl;
    calc_out << "isotope.range" << "\t" << "ion.index" << "\t" << "ion.name" << "\t" << "mz" << "\t" << "int" << "\t" << "method" << std::endl;
    scores_out << "isotope.range" << "\t" << "ion.index" << "\t" << "ion.name" << "\t" << "method" << "\t" << "label" << "\t" << "x" << "\t" << "y" << std::endl;


    const Ion precursorIon = Ion(OpenMS::AASequence::fromString("[-18.010565]ELYENKPRRPYIL"), OpenMS::Residue::Full, 3);

    //create list of b and y ions
    std::vector<Ion> ionList;
    Ion::generateFragmentIons(ionList, precursorIon.sequence, precursorIon.charge);

    double isotopeStep = OpenMS::Constants::NEUTRON_MASS_U / precursorIon.charge;

    //Loop through all spectra
    //for (int specIndex = 0; specIndex < msExperimentCentroid.getNrSpectra(); ++specIndex) {
    for (int specIndex = 0; specIndex < 10; ++specIndex) {

        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumCentroid = msExperimentCentroid.getSpectrum(specIndex);
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumProfile = msExperimentProfile.getSpectrum(specIndex);

        if (currentSpectrumCentroid.getMSLevel() == 1) continue;

        currentSpectrumCentroid.sortByPosition();

        const OpenMS::Precursor precursorInfo = currentSpectrumCentroid.getPrecursors()[0];
        std::set<OpenMS::UInt> precursorIsotopes;

        //fill precursor isotopes vector
        SpectrumUtilities::whichPrecursorIsotopes(precursorIsotopes,
                               precursorInfo,
                               precursorIon,
                               0);

        OpenMS::UInt minIsotope = *std::min_element(precursorIsotopes.begin(), precursorIsotopes.end());
        OpenMS::UInt maxIsotope = *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end());
        //loop through each fragment ion
        for (int ionIndex = 0; ionIndex < ionList.size(); ++ionIndex) {
            //compute search peak matching tolerance
            double tol = OpenMS::Math::ppmToMass(20.0, ionList[ionIndex].monoMz);

            OpenMS::Int peakIndex = -1;
            for (int i = 0; i <= maxIsotope && peakIndex == -1; ++i)
            {
                //find nearest peak to ion mz within tolerance
                peakIndex = currentSpectrumCentroid.findNearest(ionList[ionIndex].monoMz + (isotopeStep * i), tol);
            }

            // peak not found
            if (peakIndex == -1) continue;

            std::vector<std::pair<double, double> > exactPrecursorDist;
            std::vector<std::pair<double, double> > approxPrecursorFromWeightDist;
            std::vector<std::pair<double, double> > approxFragmentFromWeightDist;
            std::vector<std::pair<double, double> > approxFragmentFromWeightAndSulfurDist;
            std::vector<std::pair<double, double> > approxFragmentSplineFromWeightDist;
            std::vector<std::pair<double, double> > approxFragmentSplineFromWeightAndSulfurDist;
            std::vector<std::pair<double, double> > exactConditionalFragmentDist;
            std::vector<std::pair<double, double> > observedDist;

            SpectrumUtilities::exactConditionalFragmentIsotopeDist(exactConditionalFragmentDist,
                                                                   precursorIsotopes,
                                                                   ionList[ionIndex],
                                                                   precursorIon.sequence,
                                                                   precursorIon.charge);

            SpectrumUtilities::approxFragmentFromWeightIsotopeDist(approxFragmentFromWeightDist,
                                                precursorIsotopes,
                                                ionList[ionIndex],
                                                precursorIon.sequence,
                                                precursorIon.charge);

            SpectrumUtilities::approxFragmentFromWeightAndSIsotopeDist(approxFragmentFromWeightAndSulfurDist,
                                                    precursorIsotopes,
                                                    ionList[ionIndex],
                                                    precursorIon.sequence,
                                                    precursorIon.charge);

            SpectrumUtilities::approxFragmentSplineFromWeightIsotopeDist(approxFragmentSplineFromWeightDist,
                                                                         precursorIsotopes,
                                                                         ionList[ionIndex],
                                                                         precursorIon.sequence,
                                                                         precursorIon.charge,
                                                                         isotopeDB);

            SpectrumUtilities::approxFragmentSplineFromWeightAndSIsotopeDist(approxFragmentSplineFromWeightAndSulfurDist,
                                                                             precursorIsotopes,
                                                                             ionList[ionIndex],
                                                                             precursorIon.sequence,
                                                                             precursorIon.charge,
                                                                             isotopeDB);



            //match theoretical distribution with observed peaks
            SpectrumUtilities::observedDistribution(observedDist, exactConditionalFragmentDist, currentSpectrumCentroid);
            //scale observed intensities across distribution
            SpectrumUtilities::scaleDistribution(observedDist);



            std::string isotope_range = std::to_string(minIsotope);
            if (precursorIsotopes.size() > 1) isotope_range += "-" + std::to_string(maxIsotope);

            std::string ion_type = (ionList[ionIndex].type == OpenMS::Residue::ResidueType::BIon) ? "B" : "Y";
            ion_type += std::to_string(ionList[ionIndex].sequence.size());
            for (int i = 0; i < ionList[ionIndex].charge; ++i) ion_type += "+";
            std::string ion_name = ion_type + " " + ionList[ionIndex].sequence.toUnmodifiedString();




            float maxIntensity = 0;
            for (auto itr = currentSpectrumProfile.begin(); itr != currentSpectrumProfile.end(); ++itr)
            {
                if (itr->getMZ() >= observedDist.front().first - 0.5 &&
                    itr->getMZ() <= observedDist.back().first + 1)
                {
                    maxIntensity = std::max(maxIntensity, itr->getIntensity());
                }
            }


            out << isotope_range << "\t" << ionIndex << "\t" << ion_name << "\t" << ionList[ionIndex].monoMz - 0.5 << "\t" << 0 << std::endl;
            for (auto itr = currentSpectrumProfile.begin(); itr != currentSpectrumProfile.end(); ++itr)
            {
                if (itr->getMZ() >= observedDist.front().first - 0.5 &&
                        itr->getMZ() <= observedDist.back().first + 1)
                {
                    out << isotope_range << "\t" << ionIndex << "\t" << ion_name << "\t" << itr->getMZ() << "\t" << itr->getIntensity()/maxIntensity << std::endl;
                }
            }
            out << isotope_range << "\t" << ionIndex << "\t" << ion_name << "\t" << ionList[ionIndex].monoMz + 3.3 << "\t" << 0 << std::endl;

            SpectrumUtilities::scaleDistribution(observedDist);

            //compute chi-squared with observed to Conditional
            double exactCondFragmentX2 = Stats::computeX2(observedDist, exactConditionalFragmentDist);
            double approxFragmentFromWeightX2 = Stats::computeX2(observedDist, approxFragmentFromWeightDist);
            double approxFragmentFromWeightAndSulfurX2 = Stats::computeX2(observedDist,
                                                                          approxFragmentFromWeightAndSulfurDist);
            double approxFragmentSplineFromWeightX2 = Stats::computeX2(observedDist,
                                                                          approxFragmentSplineFromWeightDist);
            double approxFragmentSplineFromWeightAndSulfurX2 = Stats::computeX2(observedDist,
                                                                          approxFragmentSplineFromWeightAndSulfurDist);
            normalizeDist(observedDist);
            normalizeDist(exactConditionalFragmentDist);
            normalizeDist(approxFragmentFromWeightDist);
            normalizeDist(approxFragmentFromWeightAndSulfurDist);
            normalizeDist(approxFragmentSplineFromWeightDist);
            normalizeDist(approxFragmentSplineFromWeightAndSulfurDist);

            outputDist(calc_out, exactConditionalFragmentDist, ion_name, ionIndex, isotope_range, "Exact");
            outputDist(calc_out, approxFragmentFromWeightDist, ion_name, ionIndex, isotope_range, "Averagine");
            outputDist(calc_out, approxFragmentFromWeightAndSulfurDist, ion_name, ionIndex, isotope_range, "Sulfur-specific Averagine");
            outputDist(calc_out, approxFragmentSplineFromWeightDist, ion_name, ionIndex, isotope_range, "Spline");
            outputDist(calc_out, approxFragmentSplineFromWeightAndSulfurDist, ion_name, ionIndex, isotope_range, "Sulfur-specific spline");

            outputScores(scores_out, ion_name, ionIndex, isotope_range, "Exact", exactCondFragmentX2, observedDist.front().first, 1.0);
            outputScores(scores_out, ion_name, ionIndex, isotope_range, "Averagine", approxFragmentFromWeightX2, observedDist.front().first, 0.8);
            outputScores(scores_out, ion_name, ionIndex, isotope_range, "Sulfur-specific Averagine", approxFragmentFromWeightAndSulfurX2, observedDist.front().first, 0.6);
            outputScores(scores_out, ion_name, ionIndex, isotope_range, "Spline", approxFragmentSplineFromWeightX2, observedDist.front().first, 0.4);
            outputScores(scores_out, ion_name, ionIndex, isotope_range, "Sulfur-specific spline", approxFragmentSplineFromWeightAndSulfurX2, observedDist.front().first, 0.2);
        }

    }

    return 0;
}