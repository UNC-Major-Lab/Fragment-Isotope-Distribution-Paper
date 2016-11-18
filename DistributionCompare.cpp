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
    std::cout << "usage: SpecOps data_directory input_mzML_spectra_file input_idXML_PSM_file offset_mz" << std::endl;
    std::cout << "  data_directory: full path to directory containing input data";
    std::cout << " and destination for output files" << std::endl;
    std::cout << "  input_mzML_spectra_file: input .mzML file contained in the data directory" << std::endl;
    std::cout << "  input_idXML_PSM_file: input .idXML file contained in the data directory" << std::endl;
    std::cout << "  offset_mz: precursor ion isolation window offset" << std::endl;
}

int main(int argc, char * argv[])
{
    //check for correct number of command line arguments
    if (argc != 5) {
        usage();
        return 0;
    }
    //report data directory
    const std::string dataDir = argv[1];
    std::cout << "Data directory: " << dataDir << std::endl;
    //report mzML file
    const std::string mzMLFileName = argv[2];
    std::cout << "Input mzML spectra file: " << mzMLFileName << std::endl;
    //report idXML file
    const std::string idXMLFileName = argv[3];
    std::cout << "Input idXML PSMs file: " << idXMLFileName << std::endl;
    //report isolation window offset
    const double offset = std::atof(argv[4]);
    std::cout << "Isolation window offset (mz): " << offset << std::endl;

    //Load input mzML file into MSExperiment
    OpenMS::MzMLFile mzMLDataFile;
    OpenMS::MSExperiment<OpenMS::Peak1D> msExperiment;
    std::cout << "Loading input mzML file " << mzMLFileName << "..." << std::endl;
    try {
        mzMLDataFile.load(dataDir + mzMLFileName, msExperiment);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 0;
    }

    //load input idXML file
    OpenMS::IdXMLFile peptideDataFile;
    std::vector<OpenMS::ProteinIdentification> protIDs;
    std::vector<OpenMS::PeptideIdentification> pepIDs;
    std::cout << "Loading input idXML file " << idXMLFileName << "..." << std::endl;
    try {
        peptideDataFile.load(dataDir + idXMLFileName, protIDs, pepIDs);
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
        distributionScoreFile.open(dataDir + scoreFileName);
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

    //Pearson CC statistics
    distributionScoreFile << "exactPrecursorCC\t";
    distributionScoreFile << "exactCondFragmentCC\t";
    distributionScoreFile << "approxPrecursorFromWeightCC\t";
    distributionScoreFile << "approxFragmentFromWeightCC\t";
    distributionScoreFile << "approxFragmentFromWeightAndSCC\t";
    //Chi-sqaured statistics
    distributionScoreFile << "exactPrecursorX2\t";
    distributionScoreFile << "exactCondFragmentX2\t";
    distributionScoreFile << "approxPrecursorFromWeightX2\t";
    distributionScoreFile << "approxFragmentFromWeightX2\t";
    distributionScoreFile << "approxFragmentFromWeightAndSX2\t";
    //total variation distance statistics
    distributionScoreFile << "exactPrecursorVD\t";
    distributionScoreFile << "exactCondFragmentVD\t";
    distributionScoreFile << "approxPrecursorFromWeightVD\t";
    distributionScoreFile << "approxFragmentFromWeightVD\t";
    distributionScoreFile << "approxFragmentFromWeightAndSVD\t";

    //Pearson CC statistics normalized
    distributionScoreFile << "exactPrecursorCC_norm\t";
    distributionScoreFile << "exactCondFragmentCC_norm\t";
    distributionScoreFile << "approxPrecursorFromWeightCC_norm\t";
    distributionScoreFile << "approxFragmentFromWeightCC_norm\t";
    distributionScoreFile << "approxFragmentFromWeightAndSCC_norm\t";
    //Chi-sqaured statistics normalized
    distributionScoreFile << "exactPrecursorX2_norm\t";
    distributionScoreFile << "exactCondFragmentX2_norm\t";
    distributionScoreFile << "approxPrecursorFromWeightX2_norm\t";
    distributionScoreFile << "approxFragmentFromWeightX2_norm\t";
    distributionScoreFile << "approxFragmentFromWeightAndSX2_norm\t";
    //total variation distance statistics normalized
    distributionScoreFile << "exactPrecursorVD_norm\t";
    distributionScoreFile << "exactCondFragmentVD_norm\t";
    distributionScoreFile << "approxPrecursorFromWeightVD_norm\t";
    distributionScoreFile << "approxFragmentFromWeightVD_norm\t";
    distributionScoreFile << "approxFragmentFromWeightAndSVD_norm\n";

    //output file for ion identification data
    const std::string ionFileName = "ions.out";
    std::ofstream ionFile;
    try {
        ionFile.open(dataDir + ionFileName);
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
                    double tol = SpectrumUtilities::ERROR_PPM * ionList[ionIndex].monoMz;

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

                        //vector for exact theoretical precursor isotope distribution <mz, probability>
                        std::vector<std::pair<double, double> > exactPrecursorDist;
                        //vector for exact conditional fragment isotope distribution <mz, probability>
                        std::vector<std::pair<double, double> > exactConditionalFragmentDist;

                        //vector for approx. precursor isotope distribution from peptide weight <mz, probability>
                        std::vector<std::pair<double, double> > approxPrecursorFromWeightDist;
                        //vector for approx. fragment isotope distribution from peptide weight <mz, probability>
                        std::vector<std::pair<double, double> > approxFragmentFromWeightDist;
                        //vector for approx. fragment isotope dist. from peptide weight and sulfurs <mz, probability>
                        std::vector<std::pair<double, double> > approxFragmentFromWeightAndSulfurDist;

                        //vector for observed isotope distribution <mz, intensity>
                        std::vector<std::pair<double, double> > observedDist;
                        //vector for precursor isotopes captured in isolation window
                        std::vector<OpenMS::UInt> precursorIsotopes;

                        //fill precursor isotopes vector
                        SpectrumUtilities::whichPrecursorIsotopes(precursorIsotopes,
                                               precursorInfo,
                                               precursorIon,
                                               offset);

                        //fill exact theoretical precursor isotope distribution vector
                        SpectrumUtilities::exactPrecursorIsotopeDist(exactPrecursorDist,
                                                  precursorIsotopes.back() + 1,
                                                  ionList[ionIndex]);
                        //fill exact conditional isotope distribution vector
                        SpectrumUtilities::exactConditionalFragmentIsotopeDist(exactConditionalFragmentDist,
                                                            precursorIsotopes,
                                                            ionList[ionIndex],
                                                            precursorIon.sequence,
                                                            precursorIon.charge);

                        //fill approx precursor isotope distribution
                        SpectrumUtilities::approxPrecursorFromWeightIsotopeDist(approxPrecursorFromWeightDist,
                                                             precursorIsotopes,
                                                             ionList[ionIndex]);
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

                        //match theoretical distribution with observed peaks
                        SpectrumUtilities::observedDistribution(observedDist, exactPrecursorDist, currentSpectrum);
                        //scale observed intensities across distribution
                        SpectrumUtilities::scaleDistribution(observedDist);

                        //compute pearsonCC with observed to exact precursor dist
                        double exactPrecursorCC = Stats::computeCC(observedDist, exactPrecursorDist);
                        //compute pearsonCC with observed to exact conditional fragment dist
                        double exactCondFragmentCC = Stats::computeCC(observedDist, exactConditionalFragmentDist);
                        //compute pearsonCC for approximate distributions
                        double approxPrecursorFromWeightCC = Stats::computeCC(observedDist, approxPrecursorFromWeightDist);
                        double approxFragmentFromWeightCC = Stats::computeCC(observedDist, approxFragmentFromWeightDist);
                        double approxFragmentFromWeightAndSulfurCC = Stats::computeCC(observedDist,
                                                                               approxFragmentFromWeightAndSulfurDist);

                        //compute chi-squared with observed to OpenMS
                        double exactPrecursorX2 = Stats::computeX2(observedDist, exactPrecursorDist);
                        //compute chi-squared with observed to Conditional
                        double exactCondFragmentX2 = Stats::computeX2(observedDist, exactConditionalFragmentDist);
                        //compute chi-squared for approximate distributions
                        double approxPrecursorFromWeightX2 = Stats::computeX2(observedDist, approxPrecursorFromWeightDist);
                        double approxFragmentFromWeightX2 = Stats::computeX2(observedDist, approxFragmentFromWeightDist);
                        double approxFragmentFromWeightAndSulfurX2 = Stats::computeX2(observedDist,
                                                                               approxFragmentFromWeightAndSulfurDist);

                        //compute total variation distance with observed to OpenMS
                        double exactPrecursorVD = Stats::computeVD(observedDist, exactPrecursorDist);
                        //compute total variation distance with observed to Conditional
                        double exactCondFragmentVD = Stats::computeVD(observedDist, exactConditionalFragmentDist);
                        //compute total variation distance for approximate distributions
                        double approxPrecursorFromWeightVD = Stats::computeVD(observedDist, approxPrecursorFromWeightDist);
                        double approxFragmentFromWeightVD = Stats::computeVD(observedDist, approxFragmentFromWeightDist);
                        double approxFragmentFromWeightAndSulfurVD = Stats::computeVD(observedDist,
                                                                               approxFragmentFromWeightAndSulfurDist);

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
                        ++numSearchedAtDepth[exactPrecursorDist.size()];
                        if (completeFlag) {
                            ++numCompleteDists[exactPrecursorDist.size()];
                        }

                        //write distribution results to file
                        distributionScoreFile << ionID << "\t";                           //ion ID
                        distributionScoreFile << SpectrumUtilities::scaledDistributionValid(observedDist) << "\t"; //valid distribution flag
                        distributionScoreFile << ionList[ionIndex].monoWeight << "\t";    //ion dist. mono weight
                        distributionScoreFile << ionList[ionIndex].charge << "\t";        //ion distribution charge
                        distributionScoreFile << exactPrecursorDist.size() << "\t";       //distribution search depth
                        distributionScoreFile << completeFlag << "\t";                    //complete dist. found
                        distributionScoreFile << completeAtDepth << "\t";                 //complete dist. up to depth

                        distributionScoreFile << precursorIsotopes.size() << "\t";
                        for (int j = 0; j < precursorIsotopes.size(); ++j) {
                            distributionScoreFile << precursorIsotopes[j] << "|";
                        }
                        distributionScoreFile << "\t";

                        distributionScoreFile << precursorIon.sequence.getFormula(OpenMS::Residue::Full, precursorIon.charge).
                                getNumberOf(ELEMENTS->getElement("Sulfur")) << "\t";
                        distributionScoreFile << ionList[ionIndex].sequence.getFormula(ionList[ionIndex].type,
                                                                                       ionList[ionIndex].charge).
                                getNumberOf(ELEMENTS->getElement("Sulfur")) << "\t";

                        //Pearson CC for exact and approximate distributions
                        distributionScoreFile << exactPrecursorCC << "\t";
                        distributionScoreFile << exactCondFragmentCC << "\t";
                        distributionScoreFile << approxPrecursorFromWeightCC << "\t";
                        distributionScoreFile << approxFragmentFromWeightCC << "\t";
                        distributionScoreFile << approxFragmentFromWeightAndSulfurCC << "\t";
                        //Chi-sqaured for exact and approximate distributions
                        distributionScoreFile << exactPrecursorX2 << "\t";
                        distributionScoreFile << exactCondFragmentX2 << "\t";
                        distributionScoreFile << approxPrecursorFromWeightX2 << "\t";
                        distributionScoreFile << approxFragmentFromWeightX2 << "\t";
                        distributionScoreFile << approxFragmentFromWeightAndSulfurX2 << "\t";
                        //Total variation distance for exact and approximate distributions
                        distributionScoreFile << exactPrecursorVD << "\t";
                        distributionScoreFile << exactCondFragmentVD << "\t";
                        distributionScoreFile << approxPrecursorFromWeightVD << "\t";
                        distributionScoreFile << approxFragmentFromWeightVD << "\t";
                        distributionScoreFile << approxFragmentFromWeightAndSulfurVD << "\t";

                        ////////////////////////////////////////////////////////////////////////
                        //Normalize predicted distributions to observed depth, recompute and output to file
                        ///////////////////////////////////////////////////////////////////////

                        //vector for exact theoretical precursor isotope distribution <mz, probability>
                        std::vector<std::pair<double, double> > exactPrecursorDist_norm;
                        //vector for exact conditional fragment isotope distribution <mz, probability>
                        std::vector<std::pair<double, double> > exactConditionalFragmentDist_norm;

                        //vector for approx. precursor isotope distribution from peptide weight <mz, probability>
                        std::vector<std::pair<double, double> > approxPrecursorFromWeightDist_norm;
                        //vector for approx. fragment isotope distribution from peptide weight <mz, probability>
                        std::vector<std::pair<double, double> > approxFragmentFromWeightDist_norm;
                        //vector for approx. fragment isotope dist. from peptide weight and sulfurs <mz, probability>
                        std::vector<std::pair<double, double> > approxFragmentFromWeightAndSulfurDist_norm;

                        //vector for observed isotope distribution <mz, intensity>
                        std::vector<std::pair<double, double> > observedDist_norm;

                        //truncate vectors to include only up to complete depth
                        for (int k = 0; k < completeAtDepth; ++k) {
                            exactPrecursorDist_norm.push_back(exactPrecursorDist[k]);
                            exactConditionalFragmentDist_norm.push_back(exactConditionalFragmentDist[k]);
                            approxPrecursorFromWeightDist_norm.push_back(approxPrecursorFromWeightDist[k]);
                            approxFragmentFromWeightDist_norm.push_back(approxFragmentFromWeightDist[k]);
                            approxFragmentFromWeightAndSulfurDist_norm.push_back(approxFragmentFromWeightAndSulfurDist[k]);
                            observedDist_norm.push_back(observedDist[k]);
                        }

                        //scale distributions to re-normalize
                        SpectrumUtilities::scaleDistribution(exactPrecursorDist_norm);
                        SpectrumUtilities::scaleDistribution(exactConditionalFragmentDist_norm);
                        SpectrumUtilities::scaleDistribution(approxPrecursorFromWeightDist_norm);
                        SpectrumUtilities::scaleDistribution(approxFragmentFromWeightDist_norm);
                        SpectrumUtilities::scaleDistribution(approxFragmentFromWeightAndSulfurDist_norm);
                        SpectrumUtilities::scaleDistribution(observedDist_norm);

                        //compute pearsonCC with observed to exact precursor dist
                        double exactPrecursorCC_norm = Stats::computeCC(observedDist_norm, exactPrecursorDist_norm);
                        //compute pearsonCC with observed to exact conditional fragment dist
                        double exactCondFragmentCC_norm = Stats::computeCC(observedDist_norm, exactConditionalFragmentDist_norm);
                        //compute pearsonCC for approximate distributions
                        double approxPrecursorFromWeightCC_norm = Stats::computeCC(observedDist_norm, approxPrecursorFromWeightDist_norm);
                        double approxFragmentFromWeightCC_norm = Stats::computeCC(observedDist_norm, approxFragmentFromWeightDist_norm);
                        double approxFragmentFromWeightAndSulfurCC_norm = Stats::computeCC(observedDist_norm,
                                                                               approxFragmentFromWeightAndSulfurDist_norm);

                        //compute chi-squared with observed to OpenMS
                        double exactPrecursorX2_norm = Stats::computeX2(observedDist_norm, exactPrecursorDist_norm);
                        //compute chi-squared with observed to Conditional
                        double exactCondFragmentX2_norm = Stats::computeX2(observedDist_norm, exactConditionalFragmentDist_norm);
                        //compute chi-squared for approximate distributions
                        double approxPrecursorFromWeightX2_norm = Stats::computeX2(observedDist_norm, approxPrecursorFromWeightDist_norm);
                        double approxFragmentFromWeightX2_norm = Stats::computeX2(observedDist_norm, approxFragmentFromWeightDist_norm);
                        double approxFragmentFromWeightAndSulfurX2_norm = Stats::computeX2(observedDist_norm,
                                                                               approxFragmentFromWeightAndSulfurDist_norm);

                        //compute total variation distance with observed to OpenMS
                        double exactPrecursorVD_norm = Stats::computeVD(observedDist_norm, exactPrecursorDist_norm);
                        //compute total variation distance with observed to Conditional
                        double exactCondFragmentVD_norm = Stats::computeVD(observedDist_norm, exactConditionalFragmentDist_norm);
                        //compute total variation distance for approximate distributions
                        double approxPrecursorFromWeightVD_norm = Stats::computeVD(observedDist_norm, approxPrecursorFromWeightDist_norm);
                        double approxFragmentFromWeightVD_norm = Stats::computeVD(observedDist_norm, approxFragmentFromWeightDist_norm);
                        double approxFragmentFromWeightAndSulfurVD_norm = Stats::computeVD(observedDist_norm,
                                                                               approxFragmentFromWeightAndSulfurDist_norm);

                        //Pearson CC for exact and approximate distributions
                        distributionScoreFile << exactPrecursorCC_norm << "\t";
                        distributionScoreFile << exactCondFragmentCC_norm << "\t";
                        distributionScoreFile << approxPrecursorFromWeightCC_norm << "\t";
                        distributionScoreFile << approxFragmentFromWeightCC_norm << "\t";
                        distributionScoreFile << approxFragmentFromWeightAndSulfurCC_norm << "\t";
                        //Chi-sqaured for exact and approximate distributions
                        distributionScoreFile << exactPrecursorX2_norm << "\t";
                        distributionScoreFile << exactCondFragmentX2_norm << "\t";
                        distributionScoreFile << approxPrecursorFromWeightX2_norm << "\t";
                        distributionScoreFile << approxFragmentFromWeightX2_norm << "\t";
                        distributionScoreFile << approxFragmentFromWeightAndSulfurX2_norm << "\t";
                        //Total variation distance for exact and approximate distributions
                        distributionScoreFile << exactPrecursorVD_norm << "\t";
                        distributionScoreFile << exactCondFragmentVD_norm << "\t";
                        distributionScoreFile << approxPrecursorFromWeightVD_norm << "\t";
                        distributionScoreFile << approxFragmentFromWeightVD_norm << "\t";
                        distributionScoreFile << approxFragmentFromWeightAndSulfurVD_norm << "\n";

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
