//

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include "Ion.h"
#include "Stats.h"
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <fstream>

//global variables
const double FDR_THRESHOLD = 0.01;          //False discovery rate threshold (%/100)
const double ERROR_PPM = 20 * 0.000001;     //acquisition mz error for peak matching
const double NEUTRON_MASS = 1.008701;       //mass of a neutron
const double ISOLATION_WINDOW_MZ = 1.6;     //the isolation window used for precursor ion collection
const OpenMS::ElementDB* ELEMENTS = OpenMS::ElementDB::getInstance();   //element database

/**
 *
 * @param precursorIsotopes
 * @param precursorInfo
 * @param precursorIon
 * @param offset
 */
void whichPrecursorIsotopes(std::vector<OpenMS::UInt> &precursorIsotopes, const OpenMS::Precursor precursorInfo,
                            const Ion precursorIon, const double offset)
{
    //isolation window lower cutoff
    double lowerCutoff = precursorInfo.getMZ() - precursorInfo.getIsolationWindowLowerOffset() + offset;
    //isolation window upper cutoff
    double upperCutoff = precursorInfo.getMZ() + precursorInfo.getIsolationWindowUpperOffset() + offset;

    //distance between isotope peaks
    double isotopeStep = NEUTRON_MASS / precursorIon.charge;
    //max number of isotopes in window
    int maxSteps = int(std::ceil((upperCutoff - lowerCutoff) / isotopeStep));

    //loop through each isotope of precursor ion
    for (int m = 0; m < maxSteps; ++m) {
        //calculate current peak mz
        double peakMZ = precursorIon.monoMz + (m * isotopeStep);
        //check if peak mz falls within window
        if ((peakMZ >= lowerCutoff) && (peakMZ <= upperCutoff)) {
            OpenMS::UInt isoPeak = m;
            precursorIsotopes.push_back(isoPeak);
        }
    }
}

/**
 * Compute the exact theoretical fragment isotopic distribution based on the precursor isotope distribution calculator.
 * @param theoDist a vector to be filled with the theoretical isotopic distribution. Composed of a vector of pairs
 * <double, double> the first being the mz of each isotope, the second the probability of seeing the peak (equivalent
 * to the peak abundance within the distribution). Vector will be cleared before being filled with distribution.
 * @param searchDepth how many isotope peaks to report in the distribution. Search depth must be greater than 0. A
 * search depth of 1 reports only the monoisotopic peak. A search depth of 2 reports m0 and m1 peaks. ect.
 * @param ion the Ion from which the monoisotopic peak will be based.
 */
void exactPrecursorIsotopeDist(std::vector<std::pair<double, double> > &theoDist,
                               const int searchDepth, const Ion &ion)
{
    //search depth 0 reports all possible!
    //based on previous checks, this shouldn't be tested
    if (searchDepth == 0) {
        return;
    }

    //clear vector for distribution
    theoDist.clear();

    //compute isotopic distribution and get vector of isotope peaks
    std::vector<std::pair<OpenMS::Size, double> > theoPeakList =
            ion.formula.getIsotopeDistribution(searchDepth).getContainer();

    //ion mz
    double ionMZ = ion.monoWeight / ion.charge;

    //loop through calculated isotopic distribution, fill with actual mz values
    for (int i = 0; i < theoPeakList.size(); ++i) {

        //compute mz of isotope peak
        double isoMZ = ionMZ + ( NEUTRON_MASS / ion.charge ) * i;

        //set theoretical distribution pair
        std::pair<double, double> theo;
        theo.first = isoMZ;
        theo.second = theoPeakList[i].second;
        theoDist.push_back(theo);
    }
}

/**
 * Compute the exact theoretical fragment isotopic distribution based on the conditional fragment isotope distribution
 * calculator.
 * @param condDist a vector to be filled with the theoretical isotopic distribution. Composed of a vector of pairs
 * <double, double> the first being the mz of each isotope, the second the probability of seeing the peak (equivalent
 * to the peak abundance within the distribution). Vector will be cleared before being filled with distribution.
 * @param precursorIsotopes a vector representation of which precurosor isotopes were isolated within the ms2
 * isolation window. A vector <0, 1, 2> would represent the m0, m1, and m2 isotopes of an isotopic
 * distribution.
 * @param ion the Ion from which the monoisotopic peak will be based.
 * @param precursorSequence the amino acid sequence of the precursor peptide that was fragmented.
 * @param precursorCharge the charge of the precursor peptide that was fragmented.
 */
void exactConditionalFragmentIsotopeDist(std::vector<std::pair<double, double> > &condDist,
                                         const std::vector<OpenMS::UInt> &precursorIsotopes,
                                         const Ion &ion,
                                         const OpenMS::AASequence &precursorSequence,
                                         const OpenMS::Int &precursorCharge)
{
    //clear vector for distribution
    condDist.clear();

    //compute conditional isotopic distribution and get vector of isotope peaks
    OpenMS::EmpiricalFormula precursorFormula = precursorSequence.getFormula(OpenMS::Residue::Full, precursorCharge);
    std::vector<std::pair<OpenMS::Size, double> > condPeakList =
            ion.formula.getConditionalFragmentIsotopeDist(precursorFormula, precursorIsotopes).getContainer();

    //ion mz
    double ionMZ = ion.monoWeight / ion.charge;

    //loop through calculated isotopic distribution, fill with actual mz values
    for (int i = 0; i < condPeakList.size(); ++i) {

        //compute mz of isotope peak
        double isoMZ = ionMZ + ( NEUTRON_MASS / ion.charge ) * i;

        //set theoretical distribution pair
        std::pair<double, double> theo;
        theo.first = isoMZ;
        theo.second = condPeakList[i].second;
        condDist.push_back(theo);
    }
}

void approxPrecursorFromWeightIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                          const std::vector<OpenMS::UInt> &precursorIsotopes,
                                          const Ion &fragmentIon)
{
    //clear vector for distribution
    approxDist.clear();

    //construct distribution of depth at the maximum precursor isotope isolated
    OpenMS::IsotopeDistribution fragmentDist(precursorIsotopes.back() + 1);

    //estimate from fragment average weight
    fragmentDist.estimateFromPeptideWeight(fragmentIon.formula.getAverageWeight());
    //re-normalize distribution
    fragmentDist.renormalize();

    //get isotope vector
    std::vector<std::pair<OpenMS::Size, double> > isotopePeaks = fragmentDist.getContainer();

    //ion mz
    double ionMZ = fragmentIon.monoWeight / fragmentIon.charge;

    //loop through calculated isotopic distribution, fill with actual mz values
    for (int i = 0; i < isotopePeaks.size(); ++i) {

        //compute mz of isotope peak
        double isoMZ = ionMZ + ( NEUTRON_MASS / fragmentIon.charge ) * i;

        //set distribution pair
        std::pair<double, double> theo;
        theo.first = isoMZ;
        theo.second = isotopePeaks[i].second;
        approxDist.push_back(theo);
    }
}

void approxFragmentFromWeightIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                         const std::vector<OpenMS::UInt> &precursorIsotopes,
                                         const Ion &fragmentIon,
                                         const OpenMS::AASequence &precursorSequence,
                                         const OpenMS::Int &precursorCharge)
{
    //clear vector for distribution
    approxDist.clear();

    //precursor average weight
    double precursorAvgWeight = precursorSequence.getAverageWeight(OpenMS::Residue::Full, precursorCharge);
    //fragment average weight
    double fragmentAvgWeight = fragmentIon.formula.getAverageWeight();

    //construct distribution
    OpenMS::IsotopeDistribution fragmentDist(precursorIsotopes.back() + 1);

    //estimate approx distribution from peptide weight
    fragmentDist.estimateForFragmentFromPeptideWeight(precursorAvgWeight, fragmentAvgWeight, precursorIsotopes);
    //re-normalize distribution
    fragmentDist.renormalize();

    //get isotope vector
    std::vector<std::pair<OpenMS::Size, double> > isotopePeaks = fragmentDist.getContainer();

    //ion mz
    double ionMZ = fragmentIon.monoWeight / fragmentIon.charge;

    //loop through calculated isotopic distribution, fill with actual mz values
    for (int i = 0; i < isotopePeaks.size(); ++i) {

        //compute mz of isotope peak
        double isoMZ = ionMZ + ( NEUTRON_MASS / fragmentIon.charge ) * i;

        //set distribution pair
        std::pair<double, double> theo;
        theo.first = isoMZ;
        theo.second = isotopePeaks[i].second;
        approxDist.push_back(theo);
    }
}

void approxFragmentFromWeightAndSIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                             const std::vector<OpenMS::UInt> &precursorIsotopes,
                                             const Ion &fragmentIon,
                                             const OpenMS::AASequence &precursorSequence,
                                             const OpenMS::Int &precursorCharge)
{
    //clear vector for distribution
    approxDist.clear();

    //precursor average weight
    double precursorAvgWeight = precursorSequence.getAverageWeight(OpenMS::Residue::Full, precursorCharge);
    //precursor number of sulfurs
    int precursorSulfurs = precursorSequence.getFormula(OpenMS::Residue::Full,
                                                        precursorCharge).getNumberOf(ELEMENTS->getElement("Sulfur"));
    //fragment average weight
    double fragmentAvgWeight = fragmentIon.formula.getAverageWeight();
    //fragment number of sulfurs
    int fragmentSulfurs = fragmentIon.formula.getNumberOf(ELEMENTS->getElement("Sulfur"));

    //construct distribution
    OpenMS::IsotopeDistribution fragmentDist(precursorIsotopes.back() + 1);

    //estimate approx distribution from peptide weight
    fragmentDist.estimateForFragmentFromPeptideWeightAndS(precursorAvgWeight, precursorSulfurs,
                                                          fragmentAvgWeight, fragmentSulfurs,
                                                          precursorIsotopes);
    //re-normalize distribution
    fragmentDist.renormalize();

    //get isotope vector
    std::vector<std::pair<OpenMS::Size, double> > isotopePeaks = fragmentDist.getContainer();

    //ion mz
    double ionMZ = fragmentIon.monoWeight / fragmentIon.charge;

    //loop through calculated isotopic distribution, fill with actual mz values
    for (int i = 0; i < isotopePeaks.size(); ++i) {

        //compute mz of isotope peak
        double isoMZ = ionMZ + ( NEUTRON_MASS / fragmentIon.charge ) * i;

        //set distribution pair
        std::pair<double, double> theo;
        theo.first = isoMZ;
        theo.second = isotopePeaks[i].second;
        approxDist.push_back(theo);
    }
}

/**
 * Identifies an isotope distribution within a mass spectrum based on the theoretical distribution mz values.
 * @param obsDist a vector to be filled with the observed isotope distribution. Composed of a vector of pairs
 * <double, double> the first being the mz of each isotope, the second being the intensity of the peak. An
 * intensity of 0 means the peak was not found in the spectrum. Vector will be cleared before being filled
 * with the distribution.
 * @param theoDist the theoretical isotopic distribution of which peaks will be searched.
 * @param spec the MS2 spectrum from which peaks will be located.
 */
void observedDistribution(std::vector<std::pair<double, double> > &obsDist,
                          std::vector<std::pair<double, double> > &theoDist,
                          const OpenMS::MSSpectrum<OpenMS::Peak1D> &spec)
{
    //loop through each theoretical peak in isotopic distribution
    for (int i = 0; i < theoDist.size(); ++i) {

        //calculate search tolerance
        double tol = ERROR_PPM * theoDist[i].first;

        //find index of actual peak in spectrum
        OpenMS::Int isoPeakIndex = spec.findNearest(theoDist[i].first, tol);

        //observed isotope distribution pair
        std::pair<double, double> obs;

        if (isoPeakIndex == -1) {
            //peak not found
            obs.first = theoDist[i].first;
            obs.second = 0;
        } else {
            //peak found
            obs.first = spec[isoPeakIndex].getMZ();
            obs.second = spec[isoPeakIndex].getIntensity();
        }
        //add observed peak pair to vector of observed isotope distribution
        obsDist.push_back(obs);
    }
}

/**
 * Scales an isotopic distribution of peaks based on raw intensity to relative intensity which sum to 1 accross
 * all peaks in the distribution.
 * @param obsDist a vector of observed peaks within an isotopic distribution. Composed of a vector of pairs
 * <double, double> the first being the mz of each isotope, the second being the raw intensity of the peak.
 * Vector will be modified to contain scaled intensity values instead of raw intensity values.
 */
void scaleDistribution(std::vector<std::pair<double, double> > &obsDist)
{
    //sum intensities across distribution
    double totalIntensity = 0;
    for (int i = 0; i < obsDist.size(); ++i) {
        totalIntensity += obsDist[i].second;
    }

    //compute scaled intensity and replace value
    for (int j = 0; j < obsDist.size(); ++j) {
        obsDist[j].second = obsDist[j].second / totalIntensity;
    }
}

/**
 * Check if a scaled isotope distribution is following a characteristic isotope distribution.
 * @param dist scaled distribution where peak intensities sum to 1
 * @return true if distribution follows a characterist rise/fall
 */
bool scaledDistributionValid(const std::vector<std::pair<double, double> > &dist)
{
    //distribution values decreasing flag
    bool valuesDecreasing = false;

    //loop from beginning of distribution, checking for decreasing values
    for (int i = 0; i < dist.size(); ++i) {
        //if values are increasing
        if (!valuesDecreasing) {
            //check for next peak decreasing and difference is greater than 5%
            if ( (dist[i+1].second < dist[i].second) && ((dist[i].second - dist[i+1].second) > 0.05) ) {
                //next peak is greater than 5% less than current peak, distribution is decreasing
                valuesDecreasing = true;
            }
        } else {    //values are decreasing
            //check for next peak increasing and difference is greater than 5%
            if ( (dist[i+1].second > dist[i].second) && ((dist[i+1].second - dist[i].second) > 0.05) ) {
                //distribution falling but next peak increases by more than 5%
                return false;
            }
        }
    }
    //distribution follows normal rise and fall.
    return true;
}

double computeCC(const std::vector<std::pair<double, double> > &obsDist,
                 const std::vector<std::pair<double, double> > &theoDist)
{
    //vector to hold observed proportions
    std::vector<double> obsProp;
    //vector to hold theoretical proportions
    std::vector<double> theoProp;

    //check they are both the same size
    if (obsDist.size() != theoDist.size()) {
        return 0;
    }

    //fill proportions vectors from distribution parameters
    for (int i = 0; i < obsDist.size(); ++i) {
        obsProp.push_back(obsDist[i].second);
        theoProp.push_back(theoDist[i].second);
    }

    //compute pearsons correlation coefficient
    return OpenMS::Math::pearsonCorrelationCoefficient(obsProp.begin(), obsProp.end(),
                                                       theoProp.begin(), theoProp.end());
}

double computeX2(const std::vector<std::pair<double, double> > &obsDist,
                 const std::vector<std::pair<double, double> > &theoDist)
{
    //vector to hold observed proportions
    std::vector<double> obsProp;
    //vector to hold theoretical proportions
    std::vector<double> theoProp;

    //check they are both the same size
    if (obsDist.size() != theoDist.size()) {
        return -1;
    }

    //fill proportions vectors from distribution parameters
    for (int i = 0; i < obsDist.size(); ++i) {
        obsProp.push_back(obsDist[i].second);
        theoProp.push_back(theoDist[i].second);
    }

    //compute chi squared statistic
    return Stats::chiSquared(obsProp.begin(), obsProp.end(),
                             theoProp.begin(), theoProp.end());
}

double computeVD(const std::vector<std::pair<double, double> > &obsDist,
                 const std::vector<std::pair<double, double> > &theoDist)
{
    //vector to hold observed proportions
    std::vector<double> obsProp;
    //vector to hold theoretical proportions
    std::vector<double> theoProp;

    //check they are both the same size
    if (obsDist.size() != theoDist.size()) {
        return -1;
    }

    //fill proportions vectors from distribution parameters
    for (int i = 0; i < obsDist.size(); ++i) {
        obsProp.push_back(obsDist[i].second);
        theoProp.push_back(theoDist[i].second);
    }

    //compute total variation distance
    return Stats::totalVariationDistance(obsProp.begin(), obsProp.end(),
                                         theoProp.begin(), theoProp.end());
}

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
                if (std::abs(precursorInfo.getMZ() - precursorIon.monoMz) > (precursorIon.monoMz * ERROR_PPM)) {
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
                    double tol = ERROR_PPM * ionList[ionIndex].monoMz;

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
                        whichPrecursorIsotopes(precursorIsotopes,
                                               precursorInfo,
                                               precursorIon,
                                               offset);

                        //fill exact theoretical precursor isotope distribution vector
                        exactPrecursorIsotopeDist(exactPrecursorDist,
                                                  precursorIsotopes.back() + 1,
                                                  ionList[ionIndex]);
                        //fill exact conditional isotope distribution vector
                        exactConditionalFragmentIsotopeDist(exactConditionalFragmentDist,
                                                            precursorIsotopes,
                                                            ionList[ionIndex],
                                                            precursorIon.sequence,
                                                            precursorIon.charge);

                        //fill approx precursor isotope distribution
                        approxPrecursorFromWeightIsotopeDist(approxPrecursorFromWeightDist,
                                                             precursorIsotopes,
                                                             ionList[ionIndex]);
                        //fill approx fragment isotope distribution
                        approxFragmentFromWeightIsotopeDist(approxFragmentFromWeightDist,
                                                            precursorIsotopes,
                                                            ionList[ionIndex],
                                                            precursorIon.sequence,
                                                            precursorIon.charge);
                        //fill approx fragment isotope distribution with sulfurs
                        approxFragmentFromWeightAndSIsotopeDist(approxFragmentFromWeightAndSulfurDist,
                                                                precursorIsotopes,
                                                                ionList[ionIndex],
                                                                precursorIon.sequence,
                                                                precursorIon.charge);

                        //match theoretical distribution with observed peaks
                        observedDistribution(observedDist, exactPrecursorDist, currentSpectrum);
                        //scale observed intensities across distribution
                        scaleDistribution(observedDist);

                        //compute pearsonCC with observed to exact precursor dist
                        double exactPrecursorCC = computeCC(observedDist, exactPrecursorDist);
                        //compute pearsonCC with observed to exact conditional fragment dist
                        double exactCondFragmentCC = computeCC(observedDist, exactConditionalFragmentDist);
                        //compute pearsonCC for approximate distributions
                        double approxPrecursorFromWeightCC = computeCC(observedDist, approxPrecursorFromWeightDist);
                        double approxFragmentFromWeightCC = computeCC(observedDist, approxFragmentFromWeightDist);
                        double approxFragmentFromWeightAndSulfurCC = computeCC(observedDist,
                                                                               approxFragmentFromWeightAndSulfurDist);

                        //compute chi-squared with observed to OpenMS
                        double exactPrecursorX2 = computeX2(observedDist, exactPrecursorDist);
                        //compute chi-squared with observed to Conditional
                        double exactCondFragmentX2 = computeX2(observedDist, exactConditionalFragmentDist);
                        //compute chi-squared for approximate distributions
                        double approxPrecursorFromWeightX2 = computeX2(observedDist, approxPrecursorFromWeightDist);
                        double approxFragmentFromWeightX2 = computeX2(observedDist, approxFragmentFromWeightDist);
                        double approxFragmentFromWeightAndSulfurX2 = computeX2(observedDist,
                                                                               approxFragmentFromWeightAndSulfurDist);

                        //compute total variation distance with observed to OpenMS
                        double exactPrecursorVD = computeVD(observedDist, exactPrecursorDist);
                        //compute total variation distance with observed to Conditional
                        double exactCondFragmentVD = computeVD(observedDist, exactConditionalFragmentDist);
                        //compute total variation distance for approximate distributions
                        double approxPrecursorFromWeightVD = computeVD(observedDist, approxPrecursorFromWeightDist);
                        double approxFragmentFromWeightVD = computeVD(observedDist, approxFragmentFromWeightDist);
                        double approxFragmentFromWeightAndSulfurVD = computeVD(observedDist,
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
                        distributionScoreFile << scaledDistributionValid(observedDist) << "\t"; //valid distribution flag
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
                        scaleDistribution(exactPrecursorDist_norm);
                        scaleDistribution(exactConditionalFragmentDist_norm);
                        scaleDistribution(approxPrecursorFromWeightDist_norm);
                        scaleDistribution(approxFragmentFromWeightDist_norm);
                        scaleDistribution(approxFragmentFromWeightAndSulfurDist_norm);
                        scaleDistribution(observedDist_norm);

                        //compute pearsonCC with observed to exact precursor dist
                        double exactPrecursorCC_norm = computeCC(observedDist_norm, exactPrecursorDist_norm);
                        //compute pearsonCC with observed to exact conditional fragment dist
                        double exactCondFragmentCC_norm = computeCC(observedDist_norm, exactConditionalFragmentDist_norm);
                        //compute pearsonCC for approximate distributions
                        double approxPrecursorFromWeightCC_norm = computeCC(observedDist_norm, approxPrecursorFromWeightDist_norm);
                        double approxFragmentFromWeightCC_norm = computeCC(observedDist_norm, approxFragmentFromWeightDist_norm);
                        double approxFragmentFromWeightAndSulfurCC_norm = computeCC(observedDist_norm,
                                                                               approxFragmentFromWeightAndSulfurDist_norm);

                        //compute chi-squared with observed to OpenMS
                        double exactPrecursorX2_norm = computeX2(observedDist_norm, exactPrecursorDist_norm);
                        //compute chi-squared with observed to Conditional
                        double exactCondFragmentX2_norm = computeX2(observedDist_norm, exactConditionalFragmentDist_norm);
                        //compute chi-squared for approximate distributions
                        double approxPrecursorFromWeightX2_norm = computeX2(observedDist_norm, approxPrecursorFromWeightDist_norm);
                        double approxFragmentFromWeightX2_norm = computeX2(observedDist_norm, approxFragmentFromWeightDist_norm);
                        double approxFragmentFromWeightAndSulfurX2_norm = computeX2(observedDist_norm,
                                                                               approxFragmentFromWeightAndSulfurDist_norm);

                        //compute total variation distance with observed to OpenMS
                        double exactPrecursorVD_norm = computeVD(observedDist_norm, exactPrecursorDist_norm);
                        //compute total variation distance with observed to Conditional
                        double exactCondFragmentVD_norm = computeVD(observedDist_norm, exactConditionalFragmentDist_norm);
                        //compute total variation distance for approximate distributions
                        double approxPrecursorFromWeightVD_norm = computeVD(observedDist_norm, approxPrecursorFromWeightDist_norm);
                        double approxFragmentFromWeightVD_norm = computeVD(observedDist_norm, approxFragmentFromWeightDist_norm);
                        double approxFragmentFromWeightAndSulfurVD_norm = computeVD(observedDist_norm,
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

                        /*
                        //report complete distributions
                        if (observedDist.size() > 4 && completeFlag == 1) {

                            for (int i = 0; i < observedDist.size(); ++i) {
                                std::cout << "Obs: mz: " << observedDist[i].first;
                                std::cout << " prop: " << observedDist[i].second;
                                std::cout << " ExactPrec: mz: " << exactPrecursorDist[i].first;
                                std::cout << " prop: " << exactPrecursorDist[i].second;
                                std::cout << " ExactFrag: mz: " << exactConditionalFragmentDist[i].first;
                                std::cout << " prop: " << exactConditionalFragmentDist[i].second;
                                std::cout << " ApproxPrec: mz: " << approxPrecursorFromWeightDist[i].first;
                                std::cout << " prop: " << approxPrecursorFromWeightDist[i].second;
                                std::cout << " ApproxFrag: mz: " << approxFragmentFromWeightDist[i].first;
                                std::cout << " prop: " << approxFragmentFromWeightDist[i].second;
                                std::cout << " ApproxFragS: mz: " << approxFragmentFromWeightAndSulfurDist[i].first;
                                std::cout << " prop: " << approxFragmentFromWeightAndSulfurDist[i].second;
                                std::cout << std::endl;
                            }

                            std::cout << "obsDist size: " << observedDist.size() << std::endl;
                            std::cout << "exactPrec size: " << exactPrecursorDist.size() << std::endl;
                            std::cout << "exactFrag size: " << exactConditionalFragmentDist.size() << std::endl;
                            std::cout << "approxPrec size: " << approxPrecursorFromWeightDist.size() << std::endl;
                            std::cout << "approxFrag size: " << approxFragmentFromWeightDist.size() << std::endl;
                            std::cout << "approxFragS size: " << approxFragmentFromWeightAndSulfurDist.size() << std::endl;
                            std::cout << "precursor charge: " << pepCharge << std::endl;
                            std::cout << "precursor isotopes size: " << precursorIsotopes.size() << std::endl;
                            std::cout << "precursor isotopes: ";
                            for (int j = 0; j < precursorIsotopes.size(); ++j) {
                                std::cout << precursorIsotopes[j] << " ";
                            }

                            std::cout << std::endl;
                            std::cout << "***************************" << std::endl;


                            //report pearsonCC
                            std::cout << "exactPrecursorCC= " << exactPrecursorCC;
                            std::cout << " exactCondFragmentCC= " << exactCondFragmentCC << std::endl;
                            std::cout << "exactPrecursorX2= " << exactPrecursorX2;
                            std::cout << " exactCondFragmentX2= " << exactCondFragmentX2 << std::endl;
                            std::cout << "***************************" << std::endl;

                        }*/
                    }
                }//ion loop
            }//peptide hit loop
        }//PSM loop
    }//spectrum loop

    /*
    for (int i = 0; i < 10; ++i) {
        std::cout << "Number of precursors at charge " << i << ": ";
        std::cout << numPrecursAtCharge[i] << std::endl;
    }
    std::cout << "Number of matched monoisotopic ions: " << numMatchedIons << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << "Number of isotope spectra searched at depth " << i << ": ";
        std::cout << numSearchedAtDepth[i] << std::endl;
    }
    for (int i = 0; i < 10; ++i) {
        std::cout << "Number of ions matched at isotope " << i << ": ";
        std::cout << numMatchedAtDepth[i] << std::endl;
    }
    for (int i = 0; i < 10; ++i) {
        std::cout << "Number of complete distributions of depth " << i << ": ";
        std::cout << numCompleteDists[i] << std::endl;
    }


    //report on peptide hits
    std::cout << "Peptide hits: " << numPeptideHits << std::endl;
    std::cout << "Peptdie hits below FDR: " << numPeptideHitsBelowFDR << std::endl;
    std::cout << "Peptide hits below FDR/peptide hits: " << numPeptideHitsBelowFDR / double(numPeptideHits) << std::endl;
    */

    //close output files
    std::cout << "Distribution comparison scorefile written to: " + scoreFileName << std::endl;
    distributionScoreFile.close();
    std::cout << "Complete ion file written to: " + ionFileName << std::endl;
    ionFile.close();

    return 0;
}
