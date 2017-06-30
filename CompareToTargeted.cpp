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
#include "IsotopeDistributions.h"
#include "CalibrationModel.h"

static const OpenMS::IsotopeSplineDB* isotopeDB = OpenMS::IsotopeSplineDB::getInstance();
const OpenMS::ElementDB* ELEMENTS = OpenMS::ElementDB::getInstance();   //element database

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

void outputDist(std::ofstream &out, std::vector<std::pair<double, double> > &dist, std::string ion_name,
                std::string isotope_range, std::string name)
{
    for (int i = 0; i < dist.size(); ++i)
    {
        out << isotope_range << "\t" << ion_name << "\t" << dist[i].first << "\t"
                 << dist[i].second << "\t" << name << std::endl;
    }
}

void outputScores(std::ofstream &out, std::string ion_name, std::string isotope_range, std::string name, double score, double x, double y)
{
    out << isotope_range << "\t" << ion_name << "\t" << name << "\t" << score << "\t" << x << "\t" << y << std::endl;
}


void usage()
{
    std::cout << "Usage: " << std::endl;
}

void writeFileHeaders(std::ofstream &out, std::ofstream &calc_out, std::ofstream &scores_out,
                      std::ofstream &distributionScoreFile, std::ofstream &isotopeScoreFile,
                      std::ofstream &distributionScoreMS1File, std::ofstream &isotopeScoreMS1File)
{

    isotopeScoreMS1File << "width\t";
    isotopeScoreMS1File << "offset\t";
    isotopeScoreMS1File << "isotope\t";
    isotopeScoreMS1File << "residualExactPrecursor\t";
    isotopeScoreMS1File << "residualCalibratedExactPrecursor\n";

    distributionScoreMS1File << "width\t";
    distributionScoreMS1File << "offset\t";
    distributionScoreMS1File << "exactPrecursorX2\t";
    distributionScoreMS1File << "calibratedExactPrecursorX2\n";



    out << "isotope.range" << "\t";
    out << "ion.name" << "\t";
    out << "mz" << "\t";
    out << "int" << std::endl;

    calc_out << "isotope.range" << "\t";
    calc_out << "ion.name" << "\t";
    calc_out << "mz" << "\t";
    calc_out << "int" << "\t";
    calc_out << "method" << std::endl;

    scores_out << "isotope.range" << "\t";
    scores_out << "ion.name" << "\t";
    scores_out << "method" << "\t";
    scores_out << "label" << "\t";
    scores_out << "x" << "\t";
    scores_out << "y" << std::endl;

    isotopeScoreFile << "scanDesc\t";
    isotopeScoreFile << "searchDepth\t";
    isotopeScoreFile << "isotope\t";
    isotopeScoreFile << "monoMz\t";
    isotopeScoreFile << "monoMass\t";
    isotopeScoreFile << "precursorMass\t";
    isotopeScoreFile << "intensity\t";
    isotopeScoreFile << "residualExactFragment\t";
    isotopeScoreFile << "residualAveragineFragment\t";
    isotopeScoreFile << "residualAveragineSulfurFragment\t";
    isotopeScoreFile << "residualExactPrecursor\t";
    isotopeScoreFile << "residualSplineFragment\t";
    isotopeScoreFile << "residualCalibratedExactFragment\t";
    isotopeScoreFile << "residualSplineSulfurFragment\n";


    //headers for distribution score output file
    distributionScoreFile << "scanDesc\t";
    distributionScoreFile << "ionID\t";                       //ion ID of monoisotopic ion
    distributionScoreFile << "distributionValid\t";           //check for valid distribution
    distributionScoreFile << "precursorMonoWeight\t";
    distributionScoreFile << "distributionMonoWeight\t";      //ion distribution monoisotopic weight
    distributionScoreFile << "ionCharge\t";                   //ion distribtuion charge
    distributionScoreFile << "searchDepth\t";                 //distribution search depth
    distributionScoreFile << "completeFlag\t";                //full ion distribution identified in spectra
    distributionScoreFile << "completeAtDepth\t";             //ion distribution complete up to depth

    distributionScoreFile << "numPrecursorIsotopes\t";
    distributionScoreFile << "precursorIsotopes\t";

    distributionScoreFile << "precursorSulfurs\t";
    distributionScoreFile << "fragmentSulfurs\t";

    //Chi-squared statistics
    distributionScoreFile << "exactCondFragmentX2\t";
    distributionScoreFile << "approxFragmentFromWeightX2\t";
    distributionScoreFile << "approxFragmentFromWeightAndSX2\t";
    distributionScoreFile << "exactPrecursorX2\t";
    distributionScoreFile << "approxPrecursorX2\t";
    distributionScoreFile << "splineFragmentX2\t";
    distributionScoreFile << "calibratedExactFragmentX2\t";
    distributionScoreFile << "splineSulfurFragmentX2\n";
}


void calcDistributions(const Ion &precursorIon, OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumCentroid,
                       OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumProfile,
                       OpenMS::Precursor &precursorInfo,
                       std::ofstream &distributionScoreFile, std::ofstream &isotopeScoreFile, std::string scanDesc,
                       std::vector<Ion> &ionList, CalibrationModel &calModel)
{
    double isotopeStep = OpenMS::Constants::C13C12_MASSDIFF_U / precursorIon.charge;

    std::set<OpenMS::UInt> precursorIsotopes = SpectrumUtilities::whichPrecursorIsotopes(precursorInfo, precursorIon, 0);

    OpenMS::UInt minIsotope = *std::min_element(precursorIsotopes.begin(), precursorIsotopes.end());
    OpenMS::UInt maxIsotope = *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end());

    if (maxIsotope > 4) return;

    double width = precursorInfo.getIsolationWindowLowerOffset() + precursorInfo.getIsolationWindowUpperOffset();

    //loop through each fragment ion
    for (auto ion : ionList) {
        //compute search peak matching tolerance
        double tol = OpenMS::Math::ppmToMass(20.0, ion.monoMz);

        OpenMS::Int peakIndex = -1;
        for (int i = 0; i <= maxIsotope && peakIndex == -1; ++i)
        {
            //find nearest peak to ion mz within tolerance
            peakIndex = currentSpectrumCentroid.findNearest(ion.monoMz + (isotopeStep * i), tol);
        }

        // peak not found
        if (peakIndex == -1) continue;

        IsotopeDistributions isotopeDistributions(precursorIsotopes, ion, precursorIon, isotopeDB, currentSpectrumCentroid, calModel, precursorInfo, width);


        std::string isotope_range = std::to_string(minIsotope);
        if (precursorIsotopes.size() > 1) isotope_range += "-" + std::to_string(maxIsotope);

        std::string ion_type = ion.getIonType();
        std::string ion_name = ion.getIonName();


        if (!isotopeDistributions.isValid) continue;

        if (isotopeDistributions.calibratedExactFragmentX2 - isotopeDistributions.exactCondFragmentX2 > 0.2) {
            int x = 1;
        }

        //for (int i = 0; i < observedDist.size(); ++i)
        for (int i = 0; i < isotopeDistributions.completeAtDepth; ++i)
        {
            double resExactFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.exactConditionalFragmentDist[i].second;
            double resAveragineFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentFromWeightDist[i].second;
            double resAveragineSulfurFragment =
                    isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentFromWeightAndSulfurDist[i].second;
            double resExactPrecursor = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.exactPrecursorDist[i].second;
            double resSplineFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentSplineFromWeightDist[i].second;
            double resSplineSulfurFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentSplineFromWeightAndSulfurDist[i].second;
            double resCalibratedExactFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.calibratedExactConditionalFragmentDist[i].second;


            isotopeScoreFile << scanDesc << "\t" << isotopeDistributions.completeAtDepth << "\t" << i << "\t"
                             << ion.monoMz << "\t" << isotopeDistributions.observedDist[i].first << "\t"
                             << precursorIon.monoWeight << "\t"
                             << isotopeDistributions.observedDist[i].second << "\t"
                             << resExactFragment << "\t" << resAveragineFragment << "\t"
                             << resAveragineSulfurFragment << "\t" << resExactPrecursor << "\t"
                             << resSplineFragment << "\t" << resCalibratedExactFragment << "\t"
                             << resSplineSulfurFragment << std::endl;
        }


        //write distribution results to file
        distributionScoreFile << scanDesc << "\t";
        distributionScoreFile << ion_type << "\t";                           //ion ID
        distributionScoreFile << isotopeDistributions.isValid << "\t"; //valid distribution flag
        distributionScoreFile << precursorIon.monoWeight << "\t";
        distributionScoreFile << ion.monoWeight << "\t";    //ion dist. mono weight
        distributionScoreFile << ion.charge << "\t";        //ion distribution charge
        distributionScoreFile << isotopeDistributions.exactConditionalFragmentDist.size() << "\t";       //distribution search depth
        distributionScoreFile << isotopeDistributions.completeFlag << "\t";                    //complete dist. found
        distributionScoreFile << isotopeDistributions.completeAtDepth << "\t";                 //complete dist. up to depth

        distributionScoreFile << precursorIsotopes.size() << "\t";
        for (auto j : precursorIsotopes) {
            distributionScoreFile << j << "|";
        }
        distributionScoreFile << "\t";

        distributionScoreFile << precursorIon.sequence.getFormula(OpenMS::Residue::Full, precursorIon.charge).
                getNumberOf(ELEMENTS->getElement("Sulfur")) << "\t";
        distributionScoreFile << ion.sequence.getFormula(ion.type,ion.charge).
                getNumberOf(ELEMENTS->getElement("Sulfur")) << "\t";

        //Chi-squared for exact and approximate distributions
        distributionScoreFile << isotopeDistributions.exactCondFragmentX2 << "\t";
        distributionScoreFile << isotopeDistributions.approxFragmentFromWeightX2 << "\t";
        distributionScoreFile << isotopeDistributions.approxFragmentFromWeightAndSulfurX2 << "\t";
        distributionScoreFile << isotopeDistributions.exactPrecursorX2 << "\t";
        distributionScoreFile << isotopeDistributions.approxFragmentSplineFromWeightX2 << "\t";
        distributionScoreFile << isotopeDistributions.calibratedExactFragmentX2 << "\t";
        distributionScoreFile << isotopeDistributions.approxFragmentSplineFromWeightAndSulfurX2 << "\n";
    }
}



void calcSpectrumIonDistribution(Ion &precursorIon, OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumCentroid,
                                 OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumProfile,
                                 OpenMS::Precursor &precursorInfo,
                                 std::ofstream &exp_out, std::ofstream &theo_out,
                                 std::ofstream &scores_out, std::vector<Ion> ionsToPlot,
                                 CalibrationModel &calModel)
{
    double isotopeStep = OpenMS::Constants::C13C12_MASSDIFF_U / precursorIon.charge;

    std::set<OpenMS::UInt> precursorIsotopes = SpectrumUtilities::whichPrecursorIsotopes(precursorInfo, precursorIon, 0.0);

    OpenMS::UInt minIsotope = *std::min_element(precursorIsotopes.begin(), precursorIsotopes.end());
    OpenMS::UInt maxIsotope = *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end());

    double width = precursorInfo.getIsolationWindowLowerOffset() + precursorInfo.getIsolationWindowUpperOffset();

    for (auto ion : ionsToPlot) {
        //compute search peak matching tolerance
        double tol = OpenMS::Math::ppmToMass(20.0, ion.monoMz);

        OpenMS::Int peakIndex = -1;
        for (int i = 0; i <= maxIsotope && peakIndex == -1; ++i) {
            //find nearest peak to ion mz within tolerance
            peakIndex = currentSpectrumCentroid.findNearest(ion.monoMz + (isotopeStep * i), tol);
        }

        // peak not found
        if (peakIndex == -1) return;

        IsotopeDistributions isotopeDistributions(precursorIsotopes, ion, precursorIon, isotopeDB,
                                                  currentSpectrumCentroid, calModel, precursorInfo, width);

        std::string isotope_range = std::to_string(minIsotope);
        if (precursorIsotopes.size() > 1) isotope_range += "-" + std::to_string(maxIsotope);

        std::string ion_type = ion.getIonType();
        std::string ion_name = ion.getIonName();

        double minMz = ion.monoMz - 0.6;
        double maxMz = ion.monoMz + (3 / ion.charge) + 0.6;

        float maxIntensity = 0;
        for (auto itr = currentSpectrumProfile.begin(); itr != currentSpectrumProfile.end(); ++itr) {
            if (itr->getMZ() >= minMz && itr->getMZ() <= maxMz) {
                maxIntensity = std::max(maxIntensity, itr->getIntensity());
            }
        }


        // Output data for low-throughput visualization
        exp_out << isotope_range << "\t" << ion_name << "\t" << minMz - 0.61 << "\t" << 0 << std::endl;
        for (auto itr = currentSpectrumProfile.begin(); itr != currentSpectrumProfile.end(); ++itr) {
            if (itr->getMZ() >= minMz && itr->getMZ() <= maxMz) {
                exp_out << isotope_range << "\t" << ion_name << "\t" << itr->getMZ() << "\t"
                        << itr->getIntensity() / maxIntensity << std::endl;
            }
        }
        exp_out << isotope_range << "\t" << ion_name << "\t" << maxMz + 0.61 << "\t" << 0 << std::endl;


        normalizeDist(isotopeDistributions.observedDist);
        normalizeDist(isotopeDistributions.exactConditionalFragmentDist);
        normalizeDist(isotopeDistributions.approxFragmentFromWeightDist);
        normalizeDist(isotopeDistributions.approxFragmentFromWeightAndSulfurDist);
        normalizeDist(isotopeDistributions.approxFragmentSplineFromWeightDist);
        normalizeDist(isotopeDistributions.approxFragmentSplineFromWeightAndSulfurDist);

        outputDist(theo_out, isotopeDistributions.exactConditionalFragmentDist, ion_name, isotope_range, "Exact");
        outputDist(theo_out, isotopeDistributions.approxFragmentFromWeightDist, ion_name, isotope_range, "Averagine");
        outputDist(theo_out, isotopeDistributions.approxFragmentFromWeightAndSulfurDist, ion_name, isotope_range,
                   "Sulfur-specific Averagine");
        outputDist(theo_out, isotopeDistributions.approxFragmentSplineFromWeightDist, ion_name, isotope_range,
                   "Spline");
        outputDist(theo_out, isotopeDistributions.approxFragmentSplineFromWeightAndSulfurDist, ion_name, isotope_range,
                   "Sulfur-specific spline");

        outputScores(scores_out, ion_name, isotope_range, "Exact", isotopeDistributions.exactCondFragmentX2, minMz, 1.0);
        outputScores(scores_out, ion_name, isotope_range, "Averagine", isotopeDistributions.approxFragmentFromWeightX2, minMz, 0.8);
        outputScores(scores_out, ion_name, isotope_range, "Sulfur-specific Averagine", isotopeDistributions.approxFragmentFromWeightAndSulfurX2, minMz, 0.6);
        outputScores(scores_out, ion_name, isotope_range, "Spline", isotopeDistributions.approxFragmentSplineFromWeightX2, minMz, 0.4);
        outputScores(scores_out, ion_name, isotope_range, "Sulfur-specific spline", isotopeDistributions.approxFragmentSplineFromWeightAndSulfurX2, minMz, 0.2);
    }
}

std::set<int> getRepresentativeScanIndexes()
{
    std::set<int> scanIndexes;

    scanIndexes.insert(1452); // M0
    scanIndexes.insert(1459); // M1
    scanIndexes.insert(1465); // M2
    scanIndexes.insert(4328); // M3

    scanIndexes.insert(1556); // M0-M1
    scanIndexes.insert(1562); // M1-M2
    scanIndexes.insert(1568); // M2-M3

    scanIndexes.insert(1783); // M0-M2
    scanIndexes.insert(1789); // M1-M3

    scanIndexes.insert(2042); // M0-M3

    return scanIndexes;

}

void getCalibrationModel(OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumCentroid, OpenMS::Precursor &precursorInfo,
                         Ion &precursorIon, CalibrationModel &calModel)
{
    double tol = OpenMS::Math::ppmToMass(20.0, precursorIon.monoMz);
    int peakIndex = currentSpectrumCentroid.findNearest(precursorIon.monoMz, tol);
    double intensity = peakIndex == -1 ? 0 : currentSpectrumCentroid[peakIndex].getIntensity();

    double width = precursorInfo.getIsolationWindowLowerOffset() + precursorInfo.getIsolationWindowUpperOffset();
    double offset = precursorInfo.getMZ() - precursorIon.monoMz;

    calModel.addCalibrationPoint(width, offset, intensity);
}

void calcMS1Distributions(OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumCentroid, OpenMS::Precursor &precursorInfo,
                         Ion &precursorIon, CalibrationModel &calModel, std::ofstream &distributionScoreFile,
                          std::ofstream &isotopeScoreFile)
{
    double width = precursorInfo.getIsolationWindowLowerOffset() + precursorInfo.getIsolationWindowUpperOffset();
    double offset = precursorInfo.getMZ() - precursorIon.monoMz;



    std::set<OpenMS::UInt> precursorIsotopes = SpectrumUtilities::whichPrecursorIsotopes(precursorInfo, precursorIon, 0.0);

    if (precursorIsotopes.size() <= 1) return;

    OpenMS::UInt minIsotope = *std::min_element(precursorIsotopes.begin(), precursorIsotopes.end());

    IsotopeDistributions isotopeDistributions(precursorIsotopes, precursorIon, isotopeDB, currentSpectrumCentroid, calModel, width, precursorInfo);

    for (int i = 0; i < isotopeDistributions.exactPrecursorDist.size(); ++i)
    {
        double resExactPrecursor = (isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.exactPrecursorDist[i].second) / isotopeDistributions.scaledObservedDist[i].second;
        double resCalibratedExactPrecursor = (isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.calibratedExactPrecursorDist[i].second) / isotopeDistributions.scaledObservedDist[i].second;

        /*if (std::isnan(resExactPrecursor)) {
            int x = 1;
        }

        if (std::abs(width-2.4) < 0.02 && std::abs(offset+0.85) < 0.02) {
            int x = 1;
        }*/

        isotopeScoreFile << width << "\t" << offset << "\t" << minIsotope+i << "\t"
                         << resExactPrecursor << "\t"
                         << resCalibratedExactPrecursor << std::endl;
    }

    distributionScoreFile << width << "\t" << offset << "\t"
                          << isotopeDistributions.exactPrecursorX2 << "\t"
                          << isotopeDistributions.calibratedExactPrecursorX2 << std::endl;

}


int main(int argc, char * argv[])
{

    OpenMS::MzMLFile mzMLDataFileProfile, mzMLDataFileCentroid;
    OpenMS::MSExperiment<OpenMS::Peak1D> msExperimentProfile, msExperimentCentroid;

    mzMLDataFileProfile.load(argv[1], msExperimentProfile);
    mzMLDataFileCentroid.load(argv[2], msExperimentCentroid);
    // MS2
    std::ofstream exp_out(argv[3]);
    std::ofstream theo_out(argv[4]);
    std::ofstream scores_out(argv[5]);
    std::ofstream distributionScoreFile(argv[6]);
    std::ofstream isotopeScoreFile(argv[7]);
    // MS1
    std::ofstream distributionScoreMS1File(argv[8]);
    std::ofstream isotopeScoreMS1File(argv[9]);

    writeFileHeaders(exp_out, theo_out, scores_out, distributionScoreFile, isotopeScoreFile, distributionScoreMS1File, isotopeScoreMS1File);

    Ion precursorIon = Ion(OpenMS::AASequence::fromString("DRVYIHPFHL"), OpenMS::Residue::Full, 3);

    //create list of b and y ions
    std::vector<Ion> ionList = precursorIon.generateFragmentIons(0, 10000);

    std::vector<Ion> ionsToPlot;
    for (auto ion : ionList) {
        if (ion.getIonType() == "B5+" || ion.getIonType() == "B9++") {
            ionsToPlot.push_back(ion);
        }
    }

    std::set<int> representativeScanIndexes = getRepresentativeScanIndexes();

    CalibrationModel calModel;

    //Loop through all spectra
    for (int specIndex = 0; specIndex < msExperimentCentroid.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumCentroid = msExperimentCentroid.getSpectrum(specIndex);
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumProfile = msExperimentProfile.getSpectrum(specIndex);

        OpenMS::Precursor precursorInfo = currentSpectrumCentroid.getPrecursors()[0];

        double scanRange = currentSpectrumCentroid.getInstrumentSettings().getScanWindows()[0].end -
                currentSpectrumCentroid.getInstrumentSettings().getScanWindows()[0].begin;

        if (scanRange < 1000) {
            currentSpectrumCentroid.sortByPosition();
            getCalibrationModel(currentSpectrumCentroid, precursorInfo, precursorIon, calModel);
        }
    }

    calModel.finishCalibation();

    for (int specIndex = 0; specIndex < msExperimentCentroid.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumCentroid = msExperimentCentroid.getSpectrum(specIndex);
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumProfile = msExperimentProfile.getSpectrum(specIndex);

        OpenMS::Precursor precursorInfo = currentSpectrumCentroid.getPrecursors()[0];

        double scanRange = currentSpectrumCentroid.getInstrumentSettings().getScanWindows()[0].end -
                           currentSpectrumCentroid.getInstrumentSettings().getScanWindows()[0].begin;

        if (scanRange < 1000) {
            calcMS1Distributions(currentSpectrumCentroid, precursorInfo, precursorIon, calModel,
                                 distributionScoreMS1File, isotopeScoreMS1File);
        }
        else if (scanRange > 1000) {
            currentSpectrumCentroid.sortByPosition();

            if (representativeScanIndexes.find(specIndex) != representativeScanIndexes.end()) {

                calcSpectrumIonDistribution(precursorIon, currentSpectrumCentroid, currentSpectrumProfile,
                                            precursorInfo, exp_out, theo_out, scores_out, ionsToPlot, calModel);
            }

            calcDistributions(precursorIon, currentSpectrumCentroid, currentSpectrumProfile,
                              precursorInfo, distributionScoreFile, isotopeScoreFile, "test", ionsToPlot , calModel); //ionsToPlot
        }
    }

    exp_out.close();
    theo_out.close();
    scores_out.close();
    distributionScoreFile.close();
    isotopeScoreFile.close();
    distributionScoreMS1File.close();
    isotopeScoreMS1File.close();

    return 0;
}