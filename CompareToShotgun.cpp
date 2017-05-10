//
#include <fstream>

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

int main(int argc, char * argv[])
{
    //check for correct number of command line arguments
    if (argc != 5) {
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
    //Mapping statistics reported by .annotate to std out

    //reporting variables
    int ionID = 0;
    int numMatchedIons = 0;
    int numSearchedAtDepth[10] = {0};
    int numMatchedAtDepth[10] = {0};
    int numCompleteDists[10] = {0};
    int numPrecursAtCharge[10] = {0};
    int numPeptideHits = 0;
    int numPeptideHitsBelowFDR = 0;

    //output file for distribution comparison results
    const std::string scoreFileName = "distributionScores.out";
    std::ofstream distributionScoreFile;
    try {
        distributionScoreFile.open(outDir + "/" + scoreFileName);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 0;
    }
    //headers for distribution score output file
    distributionScoreFile << "ionID\t";                       //ion ID of monoisotopic ion
    distributionScoreFile << "distributionValid\t";           //check for valid distribution
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
    distributionScoreFile << "approxFragmentFromWeightAndSX2\n";
    distributionScoreFile << "exactPrecursorX2\t";
    distributionScoreFile << "approxPrecursorX2\n";

    //output file for ion identification data
    const std::string ionFileName = "ions.out";
    std::ofstream ionFile;
    try {
        ionFile.open(outDir + "/" + ionFileName);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 0;
    }
    //headers for ion identification file
    ionFile << "ionID\t";               //unique ion ID
    ionFile << "spectrumIndex\t";       //spectrum index
    ionFile << "PSMindex\t";            //PSM index
    ionFile << "peptideHitIndex\t";     //peptide hit index

    ionFile << "PSM_sequence\t";   //precursor peptide sequence
    ionFile << "PSM_charge\t";     //precursor peptide charge
    ionFile << "PSM_mz\t";

    ionFile << "precursor_charge\t";
    ionFile << "precursor_mz\t";
    ionFile << "precursor_window_low\t";
    ionFile << "precursor_window_high\t";

    ionFile << "ion_sequence\t";         //ion sequence
    ionFile << "ion_type\t";             //ion type
    ionFile << "ion_charge\t";           //ion charge
    ionFile << "ion_formula\t";       //ion formula
    ionFile << "ion_monoWeight\t";       //ion weight
    ionFile << "ion_mz\t";               //ion mz

    ionFile << "ionSearchTolerance\t";  //ion search tolerance
    ionFile << "ionFoundFlag\n";        //ion not found

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
            //get precursor informaiton
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
                if (std::abs(precursorInfo.getMZ() - precursorIon.monoMz) > (precursorIon.monoMz * SpectrumUtilities::ERROR_PPM)) {
                    //std::cout << "Warning: precursor target mz does not match PSM mz! Possible offset!" << std::endl;
                }

                //create list of b and y ions
                std::vector<Ion> ionList;
                Ion::generateFragmentIons(ionList, precursorIon.sequence, precursorIon.charge);

                //loop through each ion
                for (int ionIndex = 0; ionIndex < ionList.size(); ++ionIndex) {
                    //update ionID
                    ++ionID;

                    //compute search peak matching tolerance
                    double tol = OpenMS::Math::ppmToMass(SpectrumUtilities::ERROR_PPM, ionList[ionIndex].monoMz);

                    //find nearest peak to ion mz within tolerance
                    OpenMS::Int peakIndex = currentSpectrum.findNearest(ionList[ionIndex].monoMz, tol);

                    //write ion information to file
                    ionFile << ionID << "\t";           //unique ion ID
                    ionFile << specIndex << "\t";       //spectrum index
                    ionFile << pepIDIndex << "\t";      //PSM index
                    ionFile << pepHitIndex << "\t";     //peptide hit index

                    ionFile << precursorIon.sequence << "\t";          //precursor peptide sequence
                    ionFile << precursorIon.charge << "\t";       //precursor peptide charge
                    ionFile << precursorIon.monoMz << "\t";

                    ionFile << precursorInfo.getCharge() << "\t";
                    ionFile << precursorInfo.getMZ() << "\t";
                    ionFile << precursorInfo.getIsolationWindowLowerOffset() << "\t";
                    ionFile << precursorInfo.getIsolationWindowUpperOffset() << "\t";

                    ionFile << ionList[ionIndex].sequence << "\t";      //ion sequence
                    ionFile << ionList[ionIndex].type << "\t";          //ion type
                    ionFile << ionList[ionIndex].charge << "\t";        //ion charge
                    ionFile << ionList[ionIndex].formula << "\t";       //ion formula
                    ionFile << ionList[ionIndex].monoWeight << "\t";    //ion weight
                    ionFile << ionList[ionIndex].monoMz << "\t";              //ion mz

                    ionFile << tol << "\t";             //ion search tolerance

                    if (peakIndex == -1) {
                        //ion not found
                        //write ion information to file
                        ionFile << false << "\n";       //ion not found
                    } else {
                        //ion found
                        ++numMatchedIons;
                        //write ion information to file
                        ionFile << true << "\n";        //ion found

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
                                ++numMatchedAtDepth[i];
                                if (completeFlag) {
                                    ++completeAtDepth;
                                }
                            } else {
                                completeFlag = false;
                            }
                        }
                        ++numSearchedAtDepth[exactConditionalFragmentDist.size()];
                        if (completeAtDepth == 1) {
                            continue;
                        }
                        if (completeFlag) {
                            ++numCompleteDists[exactConditionalFragmentDist.size()];
                        }


                        //scale observed intensities across distribution
                        SpectrumUtilities::scaleDistribution(observedDist);




                        /*for (int i = 0; i < observedDist.size(); ++i) {
                            double residual = observedDist[i].second - exactConditionalFragmentDist[i].second;
                            //residual = (residual * residual) / exactConditionalFragmentDist[i].second;

                            distributionScoreFile << residual << "\t" << i << "\t" << *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) << std::endl;
                        }*/
                        //continue;

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

                        //compute chi-squared with observed to Conditional
                        double exactCondFragmentX2 = Stats::computeX2(observedDist, exactConditionalFragmentDist);
                        double approxFragmentFromWeightX2 = Stats::computeX2(observedDist, approxFragmentFromWeightDist);
                        double approxFragmentFromWeightAndSulfurX2 = Stats::computeX2(observedDist,
                                                                               approxFragmentFromWeightAndSulfurDist);
                        double exactPrecursorX2 = Stats::computeX2(observedDist, exactPrecursorDist);
                        double approxPrecursorX2 = Stats::computeX2(observedDist, approxPrecursorDist);


                        //write distribution results to file
                        distributionScoreFile << ionID << "\t";                           //ion ID
                        distributionScoreFile << SpectrumUtilities::scaledDistributionValid(observedDist) << "\t"; //valid distribution flag
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
                }//ion loop
            }//peptide hit loop
        }//PSM loop
    }//spectrum loop

    //close output files
    std::cout << "Distribution comparison scorefile written to: " + scoreFileName << std::endl;
    distributionScoreFile.close();
    std::cout << "Complete ion file written to: " + ionFileName << std::endl;
    ionFile.close();

    return 0;
}
