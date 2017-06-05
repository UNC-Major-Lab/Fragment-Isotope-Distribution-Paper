//
#include <fstream>
#include <algorithm>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include "Ion.h"
#include "Stats.h"
#include "SpectrumUtilities.h"

//global variables
const double FDR_THRESHOLD = 0.01;          //False discovery rate threshold (%/100)
const OpenMS::ElementDB* ELEMENTS = OpenMS::ElementDB::getInstance();   //element database

void usage(){
    std::cout << "usage: CompareToShotgun input_mzML_spectra_file input_idXML_PSM_file offset_mz output_directory" << std::endl;
    std::cout << "\tinput_mzML_spectra_file: path to input .mzML file " << std::endl;
    std::cout << "\tinput_idXML_PSM_file: path to input .idXML file" << std::endl;
    std::cout << "\toffset_mz: precursor ion isolation window offset" << std::endl;
    std::cout << "\toutput_directory: path to output files" << std::endl;
}

std::set<Ion> getCompleteFragmentIons(const Ion &precursorIon, const OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrum,
                                                        const OpenMS::Precursor precursorInfo, double offset,
                                                        double minMz, double maxMz) {
    int ionID = 0;
    //create list of b and y ions
    std::vector<Ion> ionList;
    std::set<Ion> ionListComplete;
    Ion::generateFragmentIons(ionList, precursorIon.sequence, precursorIon.charge, minMz, maxMz);

    //loop through each ion
    for (int ionIndex = 0; ionIndex < ionList.size(); ++ionIndex) {
        ++ionID;
        //compute search peak matching tolerance
        double tol = OpenMS::Math::ppmToMass(SpectrumUtilities::ERROR_PPM, ionList[ionIndex].monoMz);

        //find nearest peak to ion mz within tolerance
        OpenMS::Int peakIndex = currentSpectrum.findNearest(ionList[ionIndex].monoMz, tol);

        if (peakIndex != -1) {
            //vector for exact conditional fragment isotope distribution <mz, probability>
            std::vector<std::pair<double, double> > exactConditionalFragmentDist;

            //vector for observed isotope distribution <mz, intensity>
            std::vector<std::pair<double, double> > observedDist;
            //vector for precursor isotopes captured in isolation window
            std::set<OpenMS::UInt> precursorIsotopes;

            //fill precursor isotopes vector
            SpectrumUtilities::whichPrecursorIsotopes(precursorIsotopes,
                                                      precursorInfo,
                                                      precursorIon,
                                                      offset);

            //fill exact conditional isotope distribution vector
            SpectrumUtilities::exactConditionalFragmentIsotopeDist(exactConditionalFragmentDist,
                                                                   precursorIsotopes,
                                                                   ionList[ionIndex],
                                                                   precursorIon.sequence,
                                                                   precursorIon.charge);


            //match theoretical distribution with observed peaks
            SpectrumUtilities::observedDistribution(observedDist, exactConditionalFragmentDist, currentSpectrum);

            //report on matched ion distribution depth
            bool completeFlag = true;
            for (int i = 0; i < observedDist.size(); ++i) {
                if (observedDist[i].second == 0) {
                    completeFlag = false;
                }
            }

            if (completeFlag) {
                ionListComplete.insert(ionList[ionIndex]);
            }

            if (*std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) > 10) {
                int x = 1;
            }

        }
    }

    return ionListComplete;
}

void calcDistributions(const Ion &precursorIon, const OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrum,
                       const OpenMS::Precursor &precursorInfo, double offset, std::ofstream &distributionScoreFile,
                       std::ofstream &isotopeScoreFile, double minMz, double maxMz, bool isSIM)
{
    int ionID = 0;
    //create list of b and y ions
    std::vector<Ion> ionList;
    Ion::generateFragmentIons(ionList, precursorIon.sequence, precursorIon.charge, minMz, maxMz);

    //loop through each ion
    for (int ionIndex = 0; ionIndex < ionList.size(); ++ionIndex) {
        ++ionID;
        //compute search peak matching tolerance
        double tol = OpenMS::Math::ppmToMass(SpectrumUtilities::ERROR_PPM, ionList[ionIndex].monoMz);

        //find nearest peak to ion mz within tolerance
        OpenMS::Int peakIndex = currentSpectrum.findNearest(ionList[ionIndex].monoMz, tol);

        if (peakIndex != -1) {
            //vector for exact conditional fragment isotope distribution <mz, probability>
            std::vector<std::pair<double, double> > exactConditionalFragmentDist;
            //vector for approx. fragment isotope distribution from peptide weight <mz, probability>
            std::vector<std::pair<double, double> > approxFragmentFromWeightDist;
            //vector for approx. fragment isotope dist. from peptide weight and sulfurs <mz, probability>
            std::vector<std::pair<double, double> > approxFragmentFromWeightAndSulfurDist;
            std::vector<std::pair<double, double> > exactPrecursorDist;
            std::vector<std::pair<double, double> > approxPrecursorDist;


            //vector for observed isotope distribution <mz, intensity>
            std::vector<std::pair<double, double> > observedDist;
            std::vector<std::pair<double, double> > oriObservedDist;
            //vector for precursor isotopes captured in isolation window
            std::set<OpenMS::UInt> precursorIsotopes;

            //fill precursor isotopes vector
            SpectrumUtilities::whichPrecursorIsotopes(precursorIsotopes,
                                                      precursorInfo,
                                                      precursorIon,
                                                      offset);

            //fill exact conditional isotope distribution vector
            SpectrumUtilities::exactConditionalFragmentIsotopeDist(exactConditionalFragmentDist,
                                                                   precursorIsotopes,
                                                                   ionList[ionIndex],
                                                                   precursorIon.sequence,
                                                                   precursorIon.charge);


            //match theoretical distribution with observed peaks
            SpectrumUtilities::observedDistribution(observedDist, exactConditionalFragmentDist, currentSpectrum);

            //report on matched ion distribution depth
            bool completeFlag = true;
            int completeAtDepth = 0;
            for (int i = 0; i < observedDist.size(); ++i) {
                if (observedDist[i].second != 0) {
                    if (completeFlag) {
                        ++completeAtDepth;
                    }
                } else {
                    completeFlag = false;
                }
            }

            //scale observed intensities across distribution
            oriObservedDist = observedDist;
            SpectrumUtilities::scaleDistribution(observedDist);

            //fill approx fragment isotope distribution
            SpectrumUtilities::approxFragmentFromWeightIsotopeDist(approxFragmentFromWeightDist,
                                                                   precursorIsotopes,
                                                                   ionList[ionIndex],
                                                                   precursorIon.sequence,
                                                                   precursorIon.charge);
            //fill approx fragment isotope distribution with sulfurs
            SpectrumUtilities::approxFragmentFromWeightAndSIsotopeDist(approxFragmentFromWeightAndSulfurDist,
                                                                       precursorIsotopes,
                                                                       ionList[ionIndex],
                                                                       precursorIon.sequence,
                                                                       precursorIon.charge);

            SpectrumUtilities::exactPrecursorIsotopeDist(exactPrecursorDist, precursorIsotopes,
                                                         ionList[ionIndex]);

            SpectrumUtilities::approxPrecursorFromWeightIsotopeDist(approxPrecursorDist, precursorIsotopes,
                                                                    ionList[ionIndex]);


            //if (completeFlag)
            if (completeAtDepth >= 2 && completeAtDepth <= 4)
            {
                for (int i = 0; i < observedDist.size(); ++i)
                {
                    double resExactFragment = observedDist[i].second - exactConditionalFragmentDist[i].second;
                    double resAveragineFragment = observedDist[i].second - approxFragmentFromWeightDist[i].second;
                    double resAveragineSulfurFragment =
                            observedDist[i].second - approxFragmentFromWeightAndSulfurDist[i].second;
                    double resExactPrecursor = observedDist[i].second - exactPrecursorDist[i].second;
                    double resAveraginePrecursor = observedDist[i].second - approxPrecursorDist[i].second;

                    isotopeScoreFile << isSIM << "\t" << completeAtDepth << "\t" << i << "\t"
                                     << ionList[ionIndex].monoMz << "\t" << observedDist[i].first << "\t"
                                     << oriObservedDist[i].second << "\t"
                                     << resExactFragment << "\t" << resAveragineFragment << "\t"
                                     << resAveragineSulfurFragment << "\t" << resExactPrecursor << "\t"
                                     << resAveraginePrecursor << std::endl;
                }
            }

            //compute chi-squared with observed to Conditional
            double exactCondFragmentX2 = Stats::computeX2(observedDist, exactConditionalFragmentDist);
            double approxFragmentFromWeightX2 = Stats::computeX2(observedDist, approxFragmentFromWeightDist);
            double approxFragmentFromWeightAndSulfurX2 = Stats::computeX2(observedDist,
                                                                          approxFragmentFromWeightAndSulfurDist);
            double exactPrecursorX2 = Stats::computeX2(observedDist, exactPrecursorDist);
            double approxPrecursorX2 = Stats::computeX2(observedDist, approxPrecursorDist);


            //write distribution results to file
            distributionScoreFile << isSIM << "\t";
            distributionScoreFile << ionID << "\t";                           //ion ID
            distributionScoreFile << SpectrumUtilities::scaledDistributionValid(observedDist) << "\t"; //valid distribution flag
            distributionScoreFile << precursorIon.monoWeight << "\t";    //ion dist. mono weight
            distributionScoreFile << ionList[ionIndex].monoWeight << "\t";    //ion dist. mono weight
            distributionScoreFile << ionList[ionIndex].charge << "\t";        //ion distribution charge
            distributionScoreFile << exactConditionalFragmentDist.size() << "\t";       //distribution search depth
            distributionScoreFile << completeFlag << "\t";                    //complete dist. found
            distributionScoreFile << completeAtDepth << "\t";                 //complete dist. up to depth

            distributionScoreFile << precursorIsotopes.size() << "\t";
            for (auto j : precursorIsotopes) {
                distributionScoreFile << j << "|";
            }
            distributionScoreFile << "\t";

            distributionScoreFile << precursorIon.sequence.getFormula(OpenMS::Residue::Full, precursorIon.charge).
                    getNumberOf(ELEMENTS->getElement("Sulfur")) << "\t";
            distributionScoreFile << ionList[ionIndex].sequence.getFormula(ionList[ionIndex].type,
                                                                           ionList[ionIndex].charge).
                    getNumberOf(ELEMENTS->getElement("Sulfur")) << "\t";

            //Chi-squared for exact and approximate distributions
            distributionScoreFile << exactCondFragmentX2 << "\t";
            distributionScoreFile << approxFragmentFromWeightX2 << "\t";
            distributionScoreFile << approxFragmentFromWeightAndSulfurX2 << "\t";
            distributionScoreFile << exactPrecursorX2 << "\t";
            distributionScoreFile << approxPrecursorX2 << "\n";
        }
    }
}

void calcDistributions(const Ion &precursorIon, const OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrum,
                       const OpenMS::Precursor &precursorInfo, double offset, std::ofstream &distributionScoreFile,
                       std::ofstream &isotopeScoreFile, bool isSIM, std::set<Ion> &ionList)
{
    int ionID = 0;
    //loop through each ion
    for (const Ion &ion : ionList) {
        ++ionID;
        //compute search peak matching tolerance
        double tol = OpenMS::Math::ppmToMass(SpectrumUtilities::ERROR_PPM, ion.monoMz);

        //find nearest peak to ion mz within tolerance
        OpenMS::Int peakIndex = currentSpectrum.findNearest(ion.monoMz, tol);

        if (peakIndex != -1) {
            //vector for exact conditional fragment isotope distribution <mz, probability>
            std::vector<std::pair<double, double> > exactConditionalFragmentDist;
            //vector for approx. fragment isotope distribution from peptide weight <mz, probability>
            std::vector<std::pair<double, double> > approxFragmentFromWeightDist;
            //vector for approx. fragment isotope dist. from peptide weight and sulfurs <mz, probability>
            std::vector<std::pair<double, double> > approxFragmentFromWeightAndSulfurDist;
            std::vector<std::pair<double, double> > exactPrecursorDist;
            std::vector<std::pair<double, double> > approxPrecursorDist;


            //vector for observed isotope distribution <mz, intensity>
            std::vector<std::pair<double, double> > oriObservedDist;
            std::vector<std::pair<double, double> > observedDist;
            //vector for precursor isotopes captured in isolation window
            std::set<OpenMS::UInt> precursorIsotopes;

            //fill precursor isotopes vector
            SpectrumUtilities::whichPrecursorIsotopes(precursorIsotopes,
                                                      precursorInfo,
                                                      precursorIon,
                                                      offset);

            //fill exact conditional isotope distribution vector
            SpectrumUtilities::exactConditionalFragmentIsotopeDist(exactConditionalFragmentDist,
                                                                   precursorIsotopes,
                                                                   ion,
                                                                   precursorIon.sequence,
                                                                   precursorIon.charge);


            //match theoretical distribution with observed peaks
            SpectrumUtilities::observedDistribution(observedDist, exactConditionalFragmentDist, currentSpectrum);

            //report on matched ion distribution depth
            bool completeFlag = true;
            int completeAtDepth = 0;
            for (int i = 0; i < observedDist.size(); ++i) {
                if (observedDist[i].second != 0) {
                    if (completeFlag) {
                        ++completeAtDepth;
                    }
                } else {
                    completeFlag = false;
                }
            }

            //scale observed intensities across distribution
            oriObservedDist = observedDist;
            SpectrumUtilities::scaleDistribution(observedDist);


            //fill approx fragment isotope distribution
            SpectrumUtilities::approxFragmentFromWeightIsotopeDist(approxFragmentFromWeightDist,
                                                                   precursorIsotopes,
                                                                   ion,
                                                                   precursorIon.sequence,
                                                                   precursorIon.charge);
            //fill approx fragment isotope distribution with sulfurs
            SpectrumUtilities::approxFragmentFromWeightAndSIsotopeDist(approxFragmentFromWeightAndSulfurDist,
                                                                       precursorIsotopes,
                                                                       ion,
                                                                       precursorIon.sequence,
                                                                       precursorIon.charge);

            SpectrumUtilities::exactPrecursorIsotopeDist(exactPrecursorDist, precursorIsotopes,
                                                         ion);

            SpectrumUtilities::approxPrecursorFromWeightIsotopeDist(approxPrecursorDist, precursorIsotopes,
                                                                    ion);

            if (*std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) > 10) {
                int x = 1;
            }
            //if (completeFlag && completeAtDepth > 1)
            if (completeAtDepth > 1)
            {
                //for (int i = 0; i < observedDist.size(); ++i)
                for (int i = 0; i < completeAtDepth; ++i)
                {
                    double resExactFragment = observedDist[i].second - exactConditionalFragmentDist[i].second;
                    double resAveragineFragment = observedDist[i].second - approxFragmentFromWeightDist[i].second;
                    double resAveragineSulfurFragment =
                            observedDist[i].second - approxFragmentFromWeightAndSulfurDist[i].second;
                    double resExactPrecursor = observedDist[i].second - exactPrecursorDist[i].second;
                    double resAveraginePrecursor = observedDist[i].second - approxPrecursorDist[i].second;

                    isotopeScoreFile << isSIM << "\t" << completeAtDepth << "\t" << i << "\t"
                                     << ion.monoMz << "\t" << observedDist[i].first << "\t"
                                     << oriObservedDist[i].second << "\t"
                                     << resExactFragment << "\t" << resAveragineFragment << "\t"
                                     << resAveragineSulfurFragment << "\t" << resExactPrecursor << "\t"
                                     << resAveraginePrecursor << std::endl;
                }
            }

            //compute chi-squared with observed to Conditional
            double exactCondFragmentX2 = Stats::computeX2(observedDist, exactConditionalFragmentDist);
            double approxFragmentFromWeightX2 = Stats::computeX2(observedDist, approxFragmentFromWeightDist);
            double approxFragmentFromWeightAndSulfurX2 = Stats::computeX2(observedDist,
                                                                          approxFragmentFromWeightAndSulfurDist);
            double exactPrecursorX2 = Stats::computeX2(observedDist, exactPrecursorDist);
            double approxPrecursorX2 = Stats::computeX2(observedDist, approxPrecursorDist);

            bool valid = SpectrumUtilities::scaledDistributionValid(observedDist);

            //write distribution results to file
            distributionScoreFile << isSIM << "\t";
            distributionScoreFile << ionID << "\t";                           //ion ID
            distributionScoreFile << valid << "\t"; //valid distribution flag
            distributionScoreFile << precursorIon.monoWeight << "\t";
            distributionScoreFile << ion.monoWeight << "\t";    //ion dist. mono weight
            distributionScoreFile << ion.charge << "\t";        //ion distribution charge
            distributionScoreFile << exactConditionalFragmentDist.size() << "\t";       //distribution search depth
            distributionScoreFile << completeFlag << "\t";                    //complete dist. found
            distributionScoreFile << completeAtDepth << "\t";                 //complete dist. up to depth

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
            distributionScoreFile << exactCondFragmentX2 << "\t";
            distributionScoreFile << approxFragmentFromWeightX2 << "\t";
            distributionScoreFile << approxFragmentFromWeightAndSulfurX2 << "\t";
            distributionScoreFile << exactPrecursorX2 << "\t";
            distributionScoreFile << approxPrecursorX2 << "\n";
        }
    }
}

void analyzeAlternatingMS2SIMExperiment(OpenMS::MzMLFile &mzMLDataFile, OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment,
                                        double offset, double minMz, double maxMz,
                                        std::ofstream &distributionScoreFile, std::ofstream &isotopeScoreFile)
{
    std::map<int, bool> scan2isSIM;
    std::map<int, OpenMS::Precursor> scan2info;
    std::set<int> matchedScans;

    std::map<int, std::pair<Ion, std::set<Ion> > > scan2ions;

    int numSIM = 0, numNotSIM = 0;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {

            currentSpectrum.sortByPosition();

            bool isSIM = currentSpectrum.getInstrumentSettings().getScanWindows()[0].begin >= minMz;
            scan2isSIM[specIndex] = isSIM;

            if (!isSIM) {

                const OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];
                const std::vector<OpenMS::PeptideIdentification> pepIDs = currentSpectrum.getPeptideIdentifications();

                scan2info[specIndex] = precursorInfo;

                for (int pepIDIndex = 0; pepIDIndex < pepIDs.size() && pepIDIndex < 1; ++pepIDIndex) {
                    for (int pepHitIndex = 0; pepHitIndex < pepIDs[pepIDIndex].getHits().size(); ++pepHitIndex) {
                        if (pepIDs[pepIDIndex].getHits()[pepHitIndex].getScore() <= FDR_THRESHOLD) {
                            Ion precursorIon = Ion(pepIDs[pepIDIndex].getHits()[pepIDIndex].getSequence(),
                                                   OpenMS::Residue::Full,
                                                   pepIDs[pepIDIndex].getHits()[pepIDIndex].getCharge());

                            std::set<Ion> frags = getCompleteFragmentIons(precursorIon, currentSpectrum, precursorInfo, offset, minMz, maxMz);

                            scan2ions[specIndex] = std::pair<Ion, std::set<Ion> >(precursorIon, frags);

                            break;
                        }
                    }
                }
            }
        }
    }


    int numfound = 0, numnotFound = 0;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {

            if (scan2isSIM[specIndex]) {

                const OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

                bool found = false;

                for (int specIndexFull = std::max(0, specIndex-10); specIndexFull < std::min(specIndex+10, (int)msExperiment.getNrSpectra()); ++specIndexFull)
                {
                    OpenMS::MSSpectrum<OpenMS::Peak1D> fullSpectrum = msExperiment.getSpectrum(specIndexFull);

                    if (matchedScans.find(specIndexFull) == matchedScans.end()
                        && scan2info.find(specIndexFull) != scan2info.end()
                        && fullSpectrum.getPrecursors()[0].getMZ() == precursorInfo.getMZ())
                    {
                        found = true;
                        matchedScans.insert(specIndexFull);

                        if (scan2ions.find(specIndexFull) != scan2ions.end()) {

                            const Ion precursorIon = scan2ions[specIndexFull].first;

                            std::set<Ion> frags2 = scan2ions[specIndexFull].second;

                            std::set<Ion> frags = getCompleteFragmentIons(precursorIon, currentSpectrum, precursorInfo,
                                                                          offset, minMz, maxMz);


                            std::set<Ion> intersection;
                            std::set_intersection(frags2.begin(), frags2.end(),
                                                  frags.begin(), frags.end(),
                                                  std::inserter(intersection, intersection.end()));

                            scan2ions[specIndexFull] = std::pair<Ion, std::set<Ion> >(precursorIon, intersection);

                            if (intersection.size() > 0) {
                                calcDistributions(precursorIon, currentSpectrum, precursorInfo, offset,
                                                  distributionScoreFile,
                                                  isotopeScoreFile, true, intersection);
                            }
                        }
                        break;
                    }
                }

                if (!found) {
                    numnotFound++;
                } else {
                    numfound++;
                }
            }
        }
    }

    std::cout << numfound << "\t" << numnotFound << std::endl;





    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {


            if (scan2ions.find(specIndex) != scan2ions.end()) {

                const OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

                bool isSIM = scan2isSIM[specIndex];

                const Ion precursorIon = scan2ions[specIndex].first;
                std::set<Ion> intersection = scan2ions[specIndex].second;

                if (intersection.size() > 0) {
                    calcDistributions(precursorIon, currentSpectrum, precursorInfo, offset, distributionScoreFile,
                                      isotopeScoreFile, isSIM, intersection);
                }
            }
        }
    }
}

void analyzeAlternatingMS2IsoExperiment(OpenMS::MzMLFile &mzMLDataFile, OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment,
                                        double offset, std::ofstream &distributionScoreFile,
                                        std::ofstream &isotopeScoreFile)
{

    std::map<int, bool> scan2isQuad;
    std::map<int, OpenMS::Precursor> scan2info;
    std::set<int> matchedScans;

    std::map<int, std::pair<Ion, std::set<Ion> > > scan2ions;

    int numQuad = 0, numIT = 0;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {

            currentSpectrum.sortByPosition();

            bool isQuad = currentSpectrum.getInstrumentSettings().getScanWindows()[0].begin >= 200;
            scan2isQuad[specIndex] = isQuad;

            if (isQuad) {
                numQuad++;

                const OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];
                const std::vector<OpenMS::PeptideIdentification> pepIDs = currentSpectrum.getPeptideIdentifications();

                scan2info[specIndex] = precursorInfo;

                for (int pepIDIndex = 0; pepIDIndex < pepIDs.size() && pepIDIndex < 1; ++pepIDIndex) {
                    for (int pepHitIndex = 0; pepHitIndex < pepIDs[pepIDIndex].getHits().size(); ++pepHitIndex) {
                        if (pepIDs[pepIDIndex].getHits()[pepHitIndex].getScore() <= FDR_THRESHOLD) {
                            Ion precursorIon = Ion(pepIDs[pepIDIndex].getHits()[pepIDIndex].getSequence(),
                                                   OpenMS::Residue::Full,
                                                   pepIDs[pepIDIndex].getHits()[pepIDIndex].getCharge());

                            std::set<Ion> frags = getCompleteFragmentIons(precursorIon, currentSpectrum, precursorInfo, offset, 0, 10000);

                            scan2ions[specIndex] = std::pair<Ion, std::set<Ion> >(precursorIon, frags);

                            break;
                        }
                    }
                }
            }
        }
    }


    int numfound = 0, numnotFound = 0;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {

            if (!scan2isQuad[specIndex]) {
                numIT++;
                const OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

                bool found = false;

                for (int specIndexFull = std::max(0, specIndex-10); specIndexFull < std::min(specIndex+10, (int)msExperiment.getNrSpectra()); ++specIndexFull)
                {
                    OpenMS::MSSpectrum<OpenMS::Peak1D> fullSpectrum = msExperiment.getSpectrum(specIndexFull);

                    if (matchedScans.find(specIndexFull) == matchedScans.end()
                        && scan2info.find(specIndexFull) != scan2info.end()
                        && fullSpectrum.getPrecursors()[0].getMZ() == precursorInfo.getMZ())
                    {
                        found = true;
                        matchedScans.insert(specIndexFull);

                        if (scan2ions.find(specIndexFull) != scan2ions.end()) {

                            const Ion precursorIon = scan2ions[specIndexFull].first;

                            std::set<Ion> frags2 = scan2ions[specIndexFull].second;

                            std::set<Ion> frags = getCompleteFragmentIons(precursorIon, currentSpectrum, precursorInfo,
                                                                          offset, 0, 10000);


                            std::set<Ion> intersection;
                            std::set_intersection(frags2.begin(), frags2.end(),
                                                  frags.begin(), frags.end(),
                                                  std::inserter(intersection, intersection.end()));

                            scan2ions[specIndexFull] = std::pair<Ion, std::set<Ion> >(precursorIon, intersection);

                            if (intersection.size() > 0) {
                                calcDistributions(precursorIon, currentSpectrum, precursorInfo, offset,
                                                  distributionScoreFile,
                                                  isotopeScoreFile, false, intersection);
                            }
                        }
                        break;
                    }
                }

                if (!found) {
                    numnotFound++;
                } else {
                    numfound++;
                }
            }
        }
    }


    std::cout << numQuad << "\t" << numIT << std::endl;
    std::cout << numfound << "\t" << numnotFound << std::endl;





    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {


            if (scan2ions.find(specIndex) != scan2ions.end()) {

                const OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

                bool isQuad = scan2isQuad[specIndex];

                const Ion precursorIon = scan2ions[specIndex].first;
                std::set<Ion> intersection = scan2ions[specIndex].second;

                if (intersection.size() > 0) {
                    calcDistributions(precursorIon, currentSpectrum, precursorInfo, offset, distributionScoreFile,
                                      isotopeScoreFile, isQuad, intersection);
                }
            }
        }
    }
}

void analyzeMS2Experiment(OpenMS::MzMLFile &mzMLDataFile, OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment,
                          double offset, std::ofstream &distributionScoreFile,
                          std::ofstream &isotopeScoreFile)
{
    //reporting variables
    int numPeptideHits = 0;
    int numPeptideHitsBelowFDR = 0;

    std::cout << "Searching for isotope distributions..." << std::endl;


    //Loop through all spectra
    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        //get copy of current spectrum
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        //sort spectrum by mz
        currentSpectrum.sortByPosition();

        //get peptide identifications
        const std::vector<OpenMS::PeptideIdentification> pepIDs = currentSpectrum.getPeptideIdentifications();

        //Loop through each peptide identification (PSM)
        for (int pepIDIndex = 0; pepIDIndex < pepIDs.size(); ++pepIDIndex) {
            //get peptide hits
            const std::vector<OpenMS::PeptideHit> pepHits = pepIDs[pepIDIndex].getHits();


            //check for more than one precursor
            if (currentSpectrum.getPrecursors().size() > 1) {
                std::cout << "Warning: more than one precursor!" << std::endl;
            }
            //get precursor information
            const OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

            //loop through each peptide hit
            for (int pepHitIndex = 0; pepHitIndex < pepHits.size(); ++pepHitIndex) {
                ++numPeptideHits;

                //if peptide score is above FDR threshold, skip to next peptide
                if (pepHits[pepHitIndex].getScore() >= FDR_THRESHOLD) {
                    continue;
                }

                ++numPeptideHitsBelowFDR;

                //get AASequence and charge from peptide hit for precursor ion
                const Ion precursorIon = Ion(pepHits[pepHitIndex].getSequence(),
                                             OpenMS::Residue::Full,
                                             pepHits[pepHitIndex].getCharge());
                //check for precursor matching PSM peptide information
                if (precursorInfo.getCharge() != precursorIon.charge) {
                    //std::cout << "Warning: precursor target charge does not match PSM charge!" << std::endl;
                }
                if (std::abs(precursorInfo.getMZ() - precursorIon.monoMz) >
                    (precursorIon.monoMz * SpectrumUtilities::ERROR_PPM)) {
                    //std::cout << "Warning: precursor target mz does not match PSM mz! Possible offset!" << std::endl;
                }

                calcDistributions(precursorIon, currentSpectrum, precursorInfo, offset, distributionScoreFile, isotopeScoreFile, 0, 3000, false);

            }//peptide hit loop
        }//PSM loop
    }//spectrum loop
}

void analyzeMS1Data(const Ion &precursorIon, const OpenMS::MSSpectrum<OpenMS::Peak1D> &ms1Spectrum)
{
   /* double tol = OpenMS::Math::ppmToMass(SpectrumUtilities::ERROR_PPM, precursorIon.monoMz);

    //find nearest peak to ion mz within tolerance
    OpenMS::Int peakIndex = ms1Spectrum.findNearest(precursorIon.monoMz, tol);

    if (peakIndex != -1) {
        std::vector<std::pair<double, double> > exactPrecursorDist;
        std::vector<std::pair<double, double> > approxPrecursorDist;
        std::vector<std::pair<double, double> > approxPrecursorAndSDist;
        std::vector<std::pair<double, double> > observedDist;

        std::set<OpenMS::UInt> precursorIsotopes;
        for (int i = 0; i < 5; ++i) precursorIsotopes.insert(i);

        SpectrumUtilities::exactPrecursorIsotopeDist(exactPrecursorDist, precursorIsotopes, precursorIon);

        SpectrumUtilities::approxPrecursorFromWeightIsotopeDist(approxPrecursorDist, precursorIsotopes, precursorIon);

        SpectrumUtilities::observedDistribution(observedDist, exactPrecursorDist, ms1Spectrum);

        SpectrumUtilities::scaleDistribution(observedDist);

        for (int i = 0; i < observedDist.size(); ++i)
        {
            double resExactPrecursor = observedDist[i].second - exactPrecursorDist[i].second;
            double resAveraginePrecursor = observedDist[i].second - approxPrecursorDist[i].second;

            isotopeScoreFile << isSIM << "\t" << completeAtDepth << "\t" << i << "\t"
                             << ionList[ionIndex].monoMz << "\t" << observedDist[i].first << "\t"
                             << oriObservedDist[i].second << "\t"
                             << resExactFragment << "\t" << resAveragineFragment << "\t"
                             << resAveragineSulfurFragment << "\t" << resExactPrecursor << "\t"
                             << resAveraginePrecursor << std::endl;
        }
    }*/
}

void writeFileHeaders(std::ofstream &distributionScoreFile, std::ofstream &isotopeScoreFile)
{
    isotopeScoreFile << "isSIM\t";
    isotopeScoreFile << "searchDepth\t";
    isotopeScoreFile << "isotope\t";
    isotopeScoreFile << "monoMz\t";
    isotopeScoreFile << "monoMass\t";
    isotopeScoreFile << "intensity\t";
    isotopeScoreFile << "residualExactFragment\t";
    isotopeScoreFile << "residualAveragineFragment\t";
    isotopeScoreFile << "residualAveragineSulfurFragment\t";
    isotopeScoreFile << "residualExactPrecursor\t";
    isotopeScoreFile << "residualAveraginePrecursor\n";


    //headers for distribution score output file
    distributionScoreFile << "isSIM\t";
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
    distributionScoreFile << "approxPrecursorX2\n";
}

int main(int argc, char * argv[])
{
    //check for correct number of command line arguments
    if (argc != 6 && argc != 8) {
        usage();
        return 0;
    }

    //report mzML file
    const std::string mzMLFilePath = argv[1];
    std::cout << "Input mzML spectra file: " << mzMLFilePath << std::endl;
    //report idXML file
    const std::string idXMLFilePath = argv[2];
    std::cout << "Input idXML PSMs file: " << idXMLFilePath << std::endl;
    //report isolation window offset
    const double offset = std::atof(argv[3]);
    std::cout << "Isolation window offset (mz): " << offset << std::endl;
    //report data directory
    const std::string outDir = argv[4];
    std::cout << "Data directory: " << outDir << std::endl;


    //output file for distribution comparison results
    const std::string scoreFileName = "distributionScores.out";
    const std::string isotopeFileName = "isotopesScores.out";
    std::ofstream distributionScoreFile, isotopeScoreFile;
    try {
        distributionScoreFile.open(outDir + "/" + scoreFileName);
        isotopeScoreFile.open(outDir + "/" + isotopeFileName);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 0;
    }

    //Load input mzML file into MSExperiment
    OpenMS::MzMLFile mzMLDataFile;
    OpenMS::MSExperiment<OpenMS::Peak1D> msExperiment;
    std::cout << "Loading input mzML file " << mzMLFilePath << "..." << std::endl;
    try {
        mzMLDataFile.load(mzMLFilePath, msExperiment);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 0;
    }

    //load input idXML file
    OpenMS::IdXMLFile peptideDataFile;
    std::vector<OpenMS::ProteinIdentification> protIDs;
    std::vector<OpenMS::PeptideIdentification> pepIDs;
    std::cout << "Loading input idXML file " << idXMLFilePath << "..." << std::endl;
    try {
        peptideDataFile.load(idXMLFilePath, protIDs, pepIDs);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 0;
    }

    //Report spectra loaded and PSMs loaded
    std::cout << "Number of spectra loaded: " << msExperiment.getNrSpectra() << std::endl;
    std::cout << "Number of peptide identifications (PSMs): " << pepIDs.size() << std::endl;

    //Map peptide identifications with spectra
    std::cout << "Mapping PSMs to associated spectra..." << std::endl;
    OpenMS::IDMapper mapper;
    try {
        mapper.annotate(msExperiment, pepIDs, protIDs);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        return 0;
    }

    std::string expType = argv[5];

    writeFileHeaders(distributionScoreFile, isotopeScoreFile);

    if (expType == "alternatingMS2SIM") {
        std::cout << "Searching for isotope distributions..." << std::endl;
        double minMz = std::atof(argv[6]);
        double maxMz = std::atof(argv[7]);
        analyzeAlternatingMS2SIMExperiment(mzMLDataFile, msExperiment, offset, minMz, maxMz, distributionScoreFile, isotopeScoreFile);
    } else if (expType == "MS2") {
        analyzeMS2Experiment(mzMLDataFile, msExperiment, offset, distributionScoreFile, isotopeScoreFile);
    } else if (expType == "alternatingMS2Iso") {
        analyzeAlternatingMS2IsoExperiment(mzMLDataFile, msExperiment, offset, distributionScoreFile, isotopeScoreFile);
    }

    //close output files
    std::cout << "Distribution comparison scorefile written to: " + scoreFileName << std::endl;
    distributionScoreFile.close();
    isotopeScoreFile.close();

    return 0;
}

