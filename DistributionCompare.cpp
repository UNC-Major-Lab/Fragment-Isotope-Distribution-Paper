//

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include "Ion.h"
#include "Stats.h"
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <fstream>

//global variables
const int maxSearchDepth = 7;       //maximum isotope distribution search depth
const double ppm = 20 * 0.000001;   //acquisition mz error for peak matching
const double neutron = 1.008701;    //mass of neutron
//const double mzWindow = 0.8;        //1/2 of the isolation window used for precursure ion

void theoreticalDistribution(std::vector<std::pair<double, double> > &theoDist, Ion &ion, const OpenMS::Int &pepCharge)
{

    //search depth of 1 reports only monoisotopic peak
    int searchDepth = pepCharge;

    //search depth 0 reports all possible!
    //based on previous checks, this shouldn't be tested
    if (searchDepth == 0) {
        return;
    }

    //compute isotopic distribution and get vector of isotope peaks
    std::vector<std::pair<OpenMS::Size, double> > theoPeakList =
            ion.formula.getIsotopeDistribution(searchDepth).getContainer();

    //ion mz
    double ionMZ = ion.monoWeight / ion.charge;

    //loop through calculated isotopic distribution, fill with actual mz values
    for (int i = 0; i < theoPeakList.size(); ++i) {

        //compute mz of isotope peak
        double isoMZ = ionMZ + ( neutron / ion.charge ) * i;

        //set theoretical distribution pair
        std::pair<double, double> theo;
        theo.first = isoMZ;
        theo.second = theoPeakList[i].second;
        theoDist.push_back(theo);
    }
}

void conditionalDistribution(std::vector<std::pair<double, double> > &condDist,
                             const Ion &ion,
                             const OpenMS::AASequence &precursorSequence,
                             const OpenMS::Int &precursorCharge)
{
    //create vector of precursor isotopes
    std::vector<OpenMS::UInt> preIsotopes;
    for (int i = 0; i < precursorCharge; ++i) {
        OpenMS::UInt isoPeak = i;
        preIsotopes.push_back(isoPeak);
    }

    //compute conditional isotopic distribution and get vector of isotope peaks
    std::vector<std::pair<OpenMS::Size, double> > condPeakList =
            ion.formula.getConditionalFragmentIsotopeDist(
                    precursorSequence.getFormula(OpenMS::Residue::Full, precursorCharge),preIsotopes).getContainer();

    //ion mz
    double ionMZ = ion.monoWeight / ion.charge;

    //loop through calculated isotopic distribution, fill with actual mz values
    for (int i = 0; i < condPeakList.size(); ++i) {

        //compute mz of isotope peak
        double isoMZ = ionMZ + ( neutron / ion.charge ) * i;

        //set theoretical distribution pair
        std::pair<double, double> theo;
        theo.first = isoMZ;
        theo.second = condPeakList[i].second;
        condDist.push_back(theo);
    }

}

void observedDistribution(std::vector<std::pair<double, double> > &obsDist,
                          std::vector<std::pair<double, double> > &theoDist,
                          const OpenMS::MSSpectrum<OpenMS::Peak1D> &spec)
{
    //loop through each theoretical peak in isotopic distribution
    for (int i = 0; i < theoDist.size(); ++i) {

        //calculate search tolerance
        double tol = ppm * theoDist[i].first;

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
        return 0;
    }

    //fill proportions vectors from distribution parameters
    for (int i = 0; i < obsDist.size(); ++i) {
        obsProp.push_back(obsDist[i].second);
        theoProp.push_back(theoDist[i].second);
    }

    //compute pearsons correlation coefficient
    return Stats::chiSquared(obsProp.begin(), obsProp.end(),
                             theoProp.begin(), theoProp.end());
}

void usage(){
    std::cout << "usage: SpecOps data_directory input_mzML_spectra_file input_idXML_PSM_file" << std::endl;
    std::cout << "  data_directory: full path to directory containing input data";
    std::cout << " and destination for output files" << std::endl;
    std::cout << "  input_mzML_spectra_file: input .mzML file contained in the data directory" << std::endl;
    std::cout << "  input_idXML_PSM_file: input .idXML file contained in the data directory" << std::endl;
}

int main(int argc, char * argv[])
{
    if (argc != 4) {
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
    int numSearchedAtDepth[maxSearchDepth] = {0};
    int numMatchedAtDepth[maxSearchDepth] = {0};
    int numCompleteDists[maxSearchDepth] = {0};
    int numPrecursAtCharge[10] = {0};

    //output file for distribution comparison results
    std::ofstream distributionScoreFile;
    try {
        distributionScoreFile.open(dataDir + "distributionScores.out");
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        usage();
        return 0;
    }
    //headers for distribution score output file
    distributionScoreFile << "ionID\t";                       //ion ID of monoisotopic ion
    distributionScoreFile << "distributionMonoWeight\t";      //ion distribution monoisotopic weight
    distributionScoreFile << "ionCharge\t";                   //ion distribtuion charge
    distributionScoreFile << "searchDepth\t";                 //distribution search depth
    distributionScoreFile << "openMSPearsonCC\t";             //pearsonCC for OpenMS distribution
    distributionScoreFile << "conditionalPearsonCC\t";        //pearsonCC for Dennis's conditional distribution
    distributionScoreFile << "openMSChiSquared\t";            //chi-squared for OpenMS distribution
    distributionScoreFile << "conditionalChiSquared\t";       //chi-squared for Dennis's cond. distribution
    distributionScoreFile << "completeFlag\t";                //full ion distribution identified in spectra
    distributionScoreFile << "completeAtDepth\n";             //ion distribution complete up to depth

    /*
    //output file for ion identification data
    std::ofstream ionFile;
    try {
        ionFile.open(dataDir + "ions.out");
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
    ionFile << "precursorSequence\t";   //precursor peptide sequence
    ionFile << "precursorCharge\t";     //precursor peptide charge
    ionFile << "ionSequence\t";         //ion sequence
    ionFile << "ionType\t";             //ion type
    ionFile << "ionCharge\t";           //ion charge
    ionFile << "ionMolFormula\t";       //ion formula
    ionFile << "ionMonoWeight\t";       //ion weight
    ionFile << "ionMZ\t";               //ion mz
    ionFile << "ionSearchTolerance\t";  //ion search tolerance
    ionFile << "ionFoundFlag\n";        //ion not found
    */

    //Loop through all spectra
    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        //sort spectrum by mz
        msExperiment.getSpectrum(specIndex).sortByPosition();

        //get copy of current spectrum
        const OpenMS::MSSpectrum<OpenMS::Peak1D> spec = msExperiment.getSpectrum(specIndex);

        //get peptide identifications
        const std::vector<OpenMS::PeptideIdentification> pepIDs = spec.getPeptideIdentifications();

        //Loop through each peptide identification (PSM)
        for (int pepIDIndex = 0; pepIDIndex < pepIDs.size(); ++pepIDIndex) {
            //get peptide hits
            const std::vector<OpenMS::PeptideHit> pepHits = pepIDs[pepIDIndex].getHits();

            //loop through each peptide hit
            for (int pepHitIndex = 0; pepHitIndex < pepHits.size(); ++pepHitIndex) {
                //get AASequence and charge from peptide hit
                const OpenMS::AASequence pepSeq = pepHits[pepHitIndex].getSequence();
                const OpenMS::Int pepCharge = pepHits[pepHitIndex].getCharge();

                //record number of precursures at each charge state
                ++numPrecursAtCharge[pepCharge];

                //if charge state 1, skip to next peptide hit
                if (pepCharge == 1) {
                    continue;
                }

                //create list of b and y ions
                std::vector<Ion> ionList;
                Ion::generateFragmentIons(ionList, pepSeq, pepCharge);

                //create array of found not found for corresponding ion list
                //std::vector<bool> ionFound;

                //loop through each ion
                for (int ionIndex = 0; ionIndex < ionList.size(); ++ionIndex) {
                    //update ionID
                    ++ionID;

                    //ion mz
                    double mz = ionList[ionIndex].monoWeight / ionList[ionIndex].charge;

                    //compute search peak matching tolerance
                    double tol = ppm * mz;

                    //find nearest peak to ion mz within tolerance
                    OpenMS::Int peakIndex = msExperiment.getSpectrum(specIndex).findNearest(mz, tol);

                    /*
                    //write ion information to file
                    ionFile << ionID << "\t";           //unique ion ID
                    ionFile << specIndex << "\t";       //spectrum index
                    ionFile << pepIDIndex << "\t";      //PSM index
                    ionFile << pepHitIndex << "\t";     //peptide hit index
                    ionFile << pepSeq << "\t";          //precursor peptide sequence
                    ionFile << pepCharge << "\t";       //precursor peptide charge
                    ionFile << ionList[ionIndex].sequence << "\t";      //ion sequence
                    ionFile << ionList[ionIndex].type << "\t";          //ion type
                    ionFile << ionList[ionIndex].charge << "\t";        //ion charge
                    ionFile << ionList[ionIndex].formula << "\t";       //ion formula
                    ionFile << ionList[ionIndex].monoWeight << "\t";    //ion weight
                    ionFile << mz << "\t";              //ion mz
                    ionFile << tol << "\t";             //ion search tolerance
                    */

                    if (peakIndex == -1) {
                        //ion not found
                        //ionFound.push_back(FALSE);
                        //write ion information to file
                        //ionFile << false << "\n";       //ion not found
                    } else {
                        //ion found
                        ++numMatchedIons;
                        //ionFound.push_back(TRUE);
                        //vector for theoretical isotope distribution <mz, probability>
                        std::vector<std::pair<double, double> > theoDist;
                        //vector for conditional fragment isotope distribution <mz, probability>
                        std::vector<std::pair<double, double> > condDist;
                        //vector for observed isotope distribution <mz, intensity>
                        std::vector<std::pair<double, double> > obsDist;

                        //fill theoretical isotope distribution vector
                        theoreticalDistribution(theoDist, ionList[ionIndex], pepCharge);

                        //fill conditional isotope distribution vector
                        conditionalDistribution(condDist, ionList[ionIndex], pepSeq, pepCharge);

                        //match theoretical distribution with observed peaks
                        observedDistribution(obsDist, theoDist, spec);

                        //scale observed intensities across distribution
                        scaleDistribution(obsDist);

                        //compute pearsonCC with observed to OpenMS
                        double openCC = computeCC(obsDist, theoDist);

                        //compute pearsonCC with observed to Conditional
                        double condCC = computeCC(obsDist, condDist);

                        //compute chi-squared with observed to OpenMS
                        double openX2 = computeX2(obsDist, theoDist);

                        //compute chi-squared with observed to Conditional
                        double condX2 = computeX2(obsDist, condDist);

                        //report on matched ion distribution depth
                        bool completeFlag = true;
                        int completeAtDepth = 0;
                        for (int i = 0; i < obsDist.size(); ++i) {
                            if (obsDist[i].second != 0) {
                                ++numMatchedAtDepth[i];
                                if (completeFlag) {
                                    ++completeAtDepth;
                                }
                            } else {
                                completeFlag = false;
                            }
                        }
                        ++numSearchedAtDepth[theoDist.size()];
                        if (completeFlag) {
                            ++numCompleteDists[theoDist.size()];
                        }

                        //write distribution results to file
                        distributionScoreFile << ionID << "\t";             //ion ID
                        distributionScoreFile << ionList[ionIndex].monoWeight << "\t";    //ion distribution monoisotopic weight
                        distributionScoreFile << ionList[ionIndex].charge << "\t";        //ion distribution charge
                        distributionScoreFile << theoDist.size() << "\t";                 //distribution search depth
                        distributionScoreFile << openCC << "\t";            //pearsonsCC for OpenMS distribution
                        distributionScoreFile << condCC << "\t";            //pearsonCC for Dennis's conditional distribution
                        distributionScoreFile << openX2 << "\t";            //chi-squared for OpenMS distribution
                        distributionScoreFile << condX2 << "\t";            //chi-squared for Dennis's cond. distribution
                        distributionScoreFile << completeFlag << "\t";      //complete distribution found
                        distributionScoreFile << completeAtDepth << "\n";   //complete distribution up to depth

                        //write ion information to file
                        //ionFile << true << "\n";        //ion found


                        //report complete distributions
                        if (completeFlag && obsDist.size() >= 4) {

                            for (int i = 0; i < theoDist.size(); ++i) {
                                std::cout << "Theo: mz: " << theoDist[i].first << " prop: " << theoDist[i].second;
                                std::cout << " Obs: mz: " << obsDist[i].first << " prop: " << obsDist[i].second;
                                std::cout << " Cond: mz: " << condDist[i].first << " prop: " << condDist[i].second;
                                std::cout << std::endl;
                            }
                            //std::cout << "OpenMS Score: " << openMSScore << std::endl;
                            //std::cout << "Cond Score: " << condScore << std::endl;

                            //report pearsonCC
                            std::cout << "openCC= " << openCC;
                            std::cout << " condCC= " << condCC << std::endl;
                            std::cout << "openX2= " << openX2;
                            std::cout << " condX2= " << condX2 << std::endl;
                            std::cout << "***************************" << std::endl;
                        }
                    }
                }//ion loop
                /*
                //report each ion found
                std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                std::cout << "Looking at spectrum: " << specIndex << std::endl;
                //report on ions found
                for (int j = 0; j < ionList.size(); ++j) {
                    //report
                    if (ionFound[j] == TRUE) {
                        std::cout << "Looking at ion# " << j;
                        std::cout << " seq: " << ionList[j].sequence;
                        std::cout << " type: " << ionList[j].type;
                        std::cout << " charge: " << ionList[j].charge;
                        std::cout << " weight: " << ionList[j].monoWeight;
                        std::cout << " mz: " << ionList[j].monoWeight / ionList[j].charge;
                        std::cout << " found: " << ionFound[j];
                        std::cout << std::endl;
                    }
                }*/
            }//peptide hit loop
        }//PSM loop
    }//spectra loop

    for (int i = 0; i < 10; ++i) {
        std::cout << "Number of precursers at charge " << i << ": ";
        std::cout << numPrecursAtCharge[i] << std::endl;
    }
    std::cout << "Number of matched monoisotopic ions: " << numMatchedIons << std::endl;
    for (int i = 0; i < maxSearchDepth; ++i) {
        std::cout << "Number of isotope spectra searched at depth " << i << ": ";
        std::cout << numSearchedAtDepth[i] << std::endl;
    }
    for (int i = 0; i < maxSearchDepth; ++i) {
        std::cout << "Number of ions matched at isotope " << i << ": ";
        std::cout << numMatchedAtDepth[i] << std::endl;
    }
    for (int i = 0; i < maxSearchDepth; ++i) {
        std::cout << "Number of complete distributions of depth " << i << ": ";
        std::cout << numCompleteDists[i] << std::endl;
    }

    //close output files
    distributionScoreFile.close();
    //ionFile.close();

    return 0;
}
