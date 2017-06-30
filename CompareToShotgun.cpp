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
#include "IsotopeDistributions.h"

//global variables
const double FDR_THRESHOLD = 0.01;          //False discovery rate threshold (%/100)
const OpenMS::ElementDB* ELEMENTS = OpenMS::ElementDB::getInstance();   //element database
static const OpenMS::IsotopeSplineDB* isotopeDB = OpenMS::IsotopeSplineDB::getInstance();

void usage()
{
    std::cout << "usage: CompareToShotgun input_mzML_spectra_file input_idXML_PSM_file offset_mz output_directory" << std::endl;
    std::cout << "\tinput_mzML_spectra_file: path to input .mzML file " << std::endl;
    std::cout << "\tinput_idXML_PSM_file: path to input .idXML file" << std::endl;
    std::cout << "\toffset_mz: precursor ion isolation window offset" << std::endl;
    std::cout << "\toutput_directory: path to output files" << std::endl;
}

std::set<Ion> getCompleteFragmentIons(Ion &precursorIon, const OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrum,
                                      const OpenMS::Precursor precursorInfo, double offset,
                                                        double minMz, double maxMz)
{
    int ionID = 0;
    //create list of b and y ions
    std::vector<Ion> ionList = precursorIon.generateFragmentIons(minMz, maxMz);
    std::set<Ion> ionListComplete;


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
            std::set<OpenMS::UInt> precursorIsotopes = SpectrumUtilities::whichPrecursorIsotopes(precursorInfo,
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

            bool completeFlag = true;
            for (int i = 0; i < observedDist.size(); ++i) {
                if (observedDist[i].second == 0) {
                    completeFlag = false;
                }
            }

            if (completeFlag) {
                ionListComplete.insert(ionList[ionIndex]);
            }
        }
    }

    return ionListComplete;
}

void calcDistributions(Ion &precursorIon, OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrum,
                       OpenMS::Precursor &precursorInfo, double offset, std::ofstream &distributionScoreFile,
                       std::ofstream &isotopeScoreFile, double minMz, double maxMz, std::string scanDesc)
{
    CalibrationModel noModel;
    int ionID = 0;
    //create list of b and y ions
    std::vector<Ion> ionList = precursorIon.generateFragmentIons(minMz, maxMz);
    double width = precursorInfo.getIsolationWindowLowerOffset() + precursorInfo.getIsolationWindowUpperOffset();

    //loop through each ion
    for (int ionIndex = 0; ionIndex < ionList.size(); ++ionIndex) {
        ++ionID;
        //compute search peak matching tolerance
        double tol = OpenMS::Math::ppmToMass(SpectrumUtilities::ERROR_PPM, ionList[ionIndex].monoMz);

        //find nearest peak to ion mz within tolerance
        OpenMS::Int peakIndex = currentSpectrum.findNearest(ionList[ionIndex].monoMz, tol);

        if (peakIndex != -1) {

            std::set<OpenMS::UInt> precursorIsotopes = SpectrumUtilities::whichPrecursorIsotopes(precursorInfo,
                                                                                                 precursorIon,
                                                                                                 offset);

            IsotopeDistributions isotopeDistributions(precursorIsotopes, ionList[ionIndex], precursorIon, isotopeDB, currentSpectrum, noModel, precursorInfo, width);

            //OpenMS::UInt max_isotope = *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end());
            //double nextMass = (ionList[ionIndex].monoWeight + (max_isotope+1)*OpenMS::Constants::C13C12_MASSDIFF_U) / ionList[ionIndex].charge;
            //peakIndex = currentSpectrum.findNearest(nextMass , tol);
            //if (peakIndex != -1) continue;

            if (!isotopeDistributions.isValid) continue;
            //if (!isotopeDistributions.completeFlag) continue;
            if (isotopeDistributions.completeAtDepth >= 2 && isotopeDistributions.completeAtDepth <= 4)
            {
                for (int i = 0; i < isotopeDistributions.scaledObservedDist.size(); ++i)
                {
                    double resExactFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.exactConditionalFragmentDist[i].second;
                    double resAveragineFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentFromWeightDist[i].second;
                    double resAveragineSulfurFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentFromWeightAndSulfurDist[i].second;
                    double resExactPrecursor = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.exactPrecursorDist[i].second;
                    double resAveraginePrecursor = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxPrecursorFromWeightDist[i].second;

                    isotopeScoreFile << scanDesc << "\t" << isotopeDistributions.completeAtDepth << "\t" << i << "\t"
                                     << ionList[ionIndex].monoMz << "\t" << isotopeDistributions.scaledObservedDist[i].first << "\t"
                                     << precursorIon.monoWeight << "\t"
                                     << isotopeDistributions.observedDist[i].second << "\t"
                                     << resExactFragment << "\t" << resAveragineFragment << "\t"
                                     << resAveragineSulfurFragment << "\t" << resExactPrecursor << "\t"
                                     << resAveraginePrecursor << std::endl;
                }
            }

            //write distribution results to file
            distributionScoreFile << scanDesc << "\t";
            distributionScoreFile << ionID << "\t";                           //ion ID
            distributionScoreFile << isotopeDistributions.isValid << "\t"; //valid distribution flag
            distributionScoreFile << precursorIon.monoWeight << "\t";    //ion dist. mono weight
            distributionScoreFile << ionList[ionIndex].monoWeight << "\t";    //ion dist. mono weight
            distributionScoreFile << ionList[ionIndex].charge << "\t";        //ion distribution charge
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
            distributionScoreFile << ionList[ionIndex].sequence.getFormula(ionList[ionIndex].type,
                                                                           ionList[ionIndex].charge).
                    getNumberOf(ELEMENTS->getElement("Sulfur")) << "\t";

            //Chi-squared for exact and approximate distributions
            distributionScoreFile << isotopeDistributions.exactCondFragmentX2 << "\t";
            distributionScoreFile << isotopeDistributions.approxFragmentFromWeightX2 << "\t";
            distributionScoreFile << isotopeDistributions.approxFragmentFromWeightAndSulfurX2 << "\t";
            distributionScoreFile << isotopeDistributions.exactPrecursorX2 << "\t";
            distributionScoreFile << isotopeDistributions.approxPrecursorX2 << "\t";
            distributionScoreFile << isotopeDistributions.approxFragmentSplineFromWeightX2 << "\t";
            distributionScoreFile << isotopeDistributions.approxFragmentSplineFromWeightAndSulfurX2 << "\n";
        }
    }
}

void calcDistributions(Ion &precursorIon, OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrum,
                       OpenMS::Precursor &precursorInfo, double offset, std::ofstream &distributionScoreFile,
                       std::ofstream &isotopeScoreFile, std::string scanDesc, std::set<Ion> &ionList)
{
    CalibrationModel noModel;
    double width = precursorInfo.getIsolationWindowLowerOffset() + precursorInfo.getIsolationWindowUpperOffset();
    int ionID = 0;
    //loop through each ion
    for (Ion ion : ionList) {
        ++ionID;
        //compute search peak matching tolerance
        double tol = OpenMS::Math::ppmToMass(SpectrumUtilities::ERROR_PPM, ion.monoMz);

        //find nearest peak to ion mz within tolerance
        OpenMS::Int peakIndex = currentSpectrum.findNearest(ion.monoMz, tol);

        if (peakIndex != -1) {
            std::set<OpenMS::UInt> precursorIsotopes = SpectrumUtilities::whichPrecursorIsotopes(precursorInfo, precursorIon, offset);

            IsotopeDistributions isotopeDistributions(precursorIsotopes, ion, precursorIon, isotopeDB, currentSpectrum, noModel, precursorInfo, width);


            //OpenMS::UInt max_isotope = *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end());
            //double nextMass = (ion.monoWeight + (max_isotope+1)*OpenMS::Constants::C13C12_MASSDIFF_U) / ion.charge;
            //peakIndex = currentSpectrum.findNearest(nextMass , tol);
            //if (peakIndex != -1) continue;


            if (!isotopeDistributions.isValid) continue;
            //if (completeFlag && completeAtDepth > 1)
            if (isotopeDistributions.completeAtDepth > 1)
            {
                //for (int i = 0; i < observedDist.size(); ++i)
                for (int i = 0; i < isotopeDistributions.completeAtDepth; ++i)
                {
                    double resExactFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.exactConditionalFragmentDist[i].second;
                    double resAveragineFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentFromWeightDist[i].second;
                    double resAveragineSulfurFragment = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxFragmentFromWeightAndSulfurDist[i].second;
                    double resExactPrecursor = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.exactPrecursorDist[i].second;
                    double resAveraginePrecursor = isotopeDistributions.scaledObservedDist[i].second - isotopeDistributions.approxPrecursorFromWeightDist[i].second;

                    isotopeScoreFile << scanDesc << "\t" << isotopeDistributions.completeAtDepth << "\t" << i << "\t"
                                     << ion.monoMz << "\t" << isotopeDistributions.scaledObservedDist[i].first << "\t"
                                     << precursorIon.monoWeight << "\t"
                                     << isotopeDistributions.observedDist[i].second << "\t"
                                     << resExactFragment << "\t" << resAveragineFragment << "\t"
                                     << resAveragineSulfurFragment << "\t" << resExactPrecursor << "\t"
                                     << resAveraginePrecursor << std::endl;
                }
            }

            //write distribution results to file
            distributionScoreFile << scanDesc << "\t";
            distributionScoreFile << ionID << "\t";                           //ion ID
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
            distributionScoreFile << isotopeDistributions.approxPrecursorX2 << "\t";
            distributionScoreFile << isotopeDistributions.approxFragmentSplineFromWeightX2 << "\t";
            distributionScoreFile << isotopeDistributions.approxFragmentSplineFromWeightAndSulfurX2 << "\n";
        }
    }
}

std::map<int, std::string> analyzeAlternatingMS2SIMExperiment(OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment,
                                                              double minMz, double maxMz)
{
    std::map<int, std::string> scan2scanDesc;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2)
        {
            bool isSIM = currentSpectrum.getInstrumentSettings().getScanWindows()[0].begin >= minMz;
            std::string scanDesc = isSIM ? "SIM " + std::to_string(minMz) + "-" + std::to_string(maxMz) : "Full";
            scan2scanDesc[specIndex] = scanDesc;
        }
    }

    return scan2scanDesc;
}


std::map<int, std::string> analyzeAlternatingMS2FragExperiment(OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment)
{
    std::map<int, std::string> scan2scanDesc;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2)
        {
            bool isCID = currentSpectrum.getPrecursors()[0].getActivationMethods().find(OpenMS::Precursor::ActivationMethod::CID) != currentSpectrum.getPrecursors()[0].getActivationMethods().end();
            std::string scanDesc = isCID ? "CID" : "HCD";
            scan2scanDesc[specIndex] = scanDesc;
        }
    }

    return scan2scanDesc;
}


std::map<int, std::string> analyzeAlternatingMS2CIDExperiment(OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment)
{
    std::map<double, int> mz2count;
    std::map<int, std::string> scan2scanDesc;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2)
        {
            double mz = currentSpectrum.getPrecursors()[0].getMZ();
            std::string scanDesc = (++mz2count[mz] % 2 == 1) ? "CID 30" : "CID 25";
            scan2scanDesc[specIndex] = scanDesc;
        }
    }

    return scan2scanDesc;
}

std::map<int, std::string> analyzeAlternatingMS2HCDExperiment(OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment)
{
    std::map<double, int> mz2count;
    std::map<int, std::string> scan2scanDesc;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2)
        {
            double mz = currentSpectrum.getPrecursors()[0].getMZ();
            std::string scanDesc = (++mz2count[mz] % 2 == 1) ? "HCD 30" : "HCD 25";
            scan2scanDesc[specIndex] = scanDesc;
        }
    }

    return scan2scanDesc;
}


std::map<int, std::string> analyzeAlternatingMS2IsoExperiment(OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment)
{
    std::map<int, std::string> scan2scanDesc;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2)
        {
            bool isQuad = currentSpectrum.getInstrumentSettings().getScanWindows()[0].begin >= 200;
            std::string scanDesc = isQuad ? "Quadrupole isolation" : "Ion Trap isolation";
            scan2scanDesc[specIndex] = scanDesc;
        }
    }

    return scan2scanDesc;
}


void analyzeMS2Experiment(OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment,
                          double offset, std::ofstream &distributionScoreFile,
                          std::ofstream &isotopeScoreFile, std::string expType)
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
            OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

            //loop through each peptide hit
            for (int pepHitIndex = 0; pepHitIndex < pepHits.size(); ++pepHitIndex) {
                ++numPeptideHits;

                //if peptide score is above FDR threshold, skip to next peptide
                if (pepHits[pepHitIndex].getScore() >= FDR_THRESHOLD) {
                    continue;
                }

                ++numPeptideHitsBelowFDR;

                //get AASequence and charge from peptide hit for precursor ion
                Ion precursorIon = Ion(pepHits[pepHitIndex].getSequence(),
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

                calcDistributions(precursorIon, currentSpectrum, precursorInfo, offset, distributionScoreFile, isotopeScoreFile, 0, 3000, expType);

            }//peptide hit loop
        }//PSM loop
    }//spectrum loop
}

void analyzeAlternatingMS2Experiment(OpenMS::MSExperiment<OpenMS::Peak1D> &msExperiment,
                                     double offset, std::ofstream &distributionScoreFile, std::ofstream &isotopeScoreFile,
                                     std::map<int, std::string> &scan2scanDesc, std::map<std::string, bool> &scanDesc2doSeq,
                                     double minMz, double maxMz)
{
    std::map<int, OpenMS::Precursor> scan2info;
    std::set<int> matchedScans;
    std::map<int, std::pair<Ion, std::set<Ion> > > scan2ions;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex)
    {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {

            currentSpectrum.sortByPosition();

            double mz = currentSpectrum.getPrecursors()[0].getMZ();
            std::string scanDesc = scan2scanDesc[specIndex];
            bool doSeq = scanDesc2doSeq[scanDesc];

            if (doSeq) {

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

            std::string scanDesc = scan2scanDesc[specIndex];
            bool doSeq = scanDesc2doSeq[scanDesc];

            if (!doSeq) {

                OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

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

                            Ion precursorIon = scan2ions[specIndexFull].first;

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
                                                  isotopeScoreFile, scanDesc, intersection);
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

    std::cout << "Found: " << numfound << std::endl;
    std::cout << "Not Found: " << numnotFound << std::endl;


    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        if (currentSpectrum.getMSLevel() == 2) {

            if (scan2ions.find(specIndex) != scan2ions.end()) {

                OpenMS::Precursor precursorInfo = currentSpectrum.getPrecursors()[0];

                std::string scanDesc = scan2scanDesc[specIndex];

                Ion precursorIon = scan2ions[specIndex].first;
                std::set<Ion> intersection = scan2ions[specIndex].second;

                if (intersection.size() > 0) {
                    calcDistributions(precursorIon, currentSpectrum, precursorInfo, offset, distributionScoreFile,
                                      isotopeScoreFile, scanDesc, intersection);
                }
            }
        }
    }
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

void countScanDesc(std::map<int, std::string> &scan2scanDesc)
{
    std::map<std::string, int> scanDesc2count;

    for (auto itr : scan2scanDesc) scanDesc2count[itr.second]++;
    for (auto itr : scanDesc2count) std::cout << itr.first << "\t" << scanDesc2count[itr.first] << std::endl;
}

void writeFileHeaders(std::ofstream &distributionScoreFile, std::ofstream &isotopeScoreFile)
{
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
    isotopeScoreFile << "residualAveraginePrecursor\n";


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
    distributionScoreFile << "approxFragmentSplineFromWeightX2\t";
    distributionScoreFile << "approxFragmentSplineFromWeightAndSulfurX2\n";
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

    std::string alternating = argv[5];
    std::string msLevel = argv[6];
    std::string expType = argv[7];

    writeFileHeaders(distributionScoreFile, isotopeScoreFile);
    std::cout << "Searching for isotope distributions..." << std::endl;

    double minMz = 0, maxMz = 10000;

    if (alternating == "alternating")
    {
        std::map<int, std::string> scan2scanDesc;
        std::map<std::string, bool> scanDesc2doSeq;

        if (msLevel == "MS2")
        {
            if (expType == "SIM_vs_Full") {
                minMz = std::atof(argv[8]);
                maxMz = std::atof(argv[9]);
                scan2scanDesc = analyzeAlternatingMS2SIMExperiment(msExperiment, minMz, maxMz);
                scanDesc2doSeq["SIM"] = false;
                scanDesc2doSeq["Full"] = true;
            } else if (expType == "Quad_vs_IT") {
                scan2scanDesc = analyzeAlternatingMS2IsoExperiment(msExperiment);
                scanDesc2doSeq["Quadrupole isolation"] = true;
                scanDesc2doSeq["Ion Trap isolation"] = false;
            } else if (expType == "HCD_vs_CID") {
                scan2scanDesc = analyzeAlternatingMS2FragExperiment(msExperiment);
                scanDesc2doSeq["CID"] = false;
                scanDesc2doSeq["HCD"] = true;
            } else if (expType == "CID30_vs_CID25") {
                scan2scanDesc = analyzeAlternatingMS2CIDExperiment(msExperiment);
                scanDesc2doSeq["CID 30"] = true;
                scanDesc2doSeq["CID 25"] = false;
            } else if (expType == "HCD30_vs_HCD25") {
                scan2scanDesc = analyzeAlternatingMS2HCDExperiment(msExperiment);
                scanDesc2doSeq["HCD 30"] = true;
                scanDesc2doSeq["HCD 25"] = false;
            }

            countScanDesc(scan2scanDesc);

            analyzeAlternatingMS2Experiment(msExperiment, offset, distributionScoreFile, isotopeScoreFile,
                                            scan2scanDesc, scanDesc2doSeq, minMz, maxMz);

        } else {

        }
    }
    else
    {
        if (msLevel == "MS2") {
            analyzeMS2Experiment(msExperiment, offset, distributionScoreFile, isotopeScoreFile, "HCD 35");
        } else {

        }
    }







    //close output files
    std::cout << "Distribution comparison scorefile written to: " + scoreFileName << std::endl;
    distributionScoreFile.close();
    isotopeScoreFile.close();

    return 0;
}

