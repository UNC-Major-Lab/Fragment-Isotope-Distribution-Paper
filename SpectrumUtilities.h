//
// Created by Dennis Goldfarb on 11/17/16.
//

#ifndef FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_SPECTRUMUTILITIES_H
#define FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_SPECTRUMUTILITIES_H

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include "Ion.h"
#include "CalibrationModel.h"



namespace SpectrumUtilities {

    static const double ERROR_PPM = 20;     //acquisition mz error for peak matching

    static const OpenMS::ElementDB* ELEMENTS = OpenMS::ElementDB::getInstance();   //element database

    std::set<OpenMS::UInt> whichPrecursorIsotopes(const OpenMS::Precursor &precursorInfo,
                                       const Ion &precursorIon, const double offset) {
        std::set<OpenMS::UInt> precursorIsotopes;
        //isolation window lower cutoff
        double lowerCutoff = precursorInfo.getMZ() - precursorInfo.getIsolationWindowLowerOffset() + offset;
        //isolation window upper cutoff
        double upperCutoff = precursorInfo.getMZ() + precursorInfo.getIsolationWindowUpperOffset() + offset;

        //distance between isotope peaks
        double isotopeStep = OpenMS::Constants::C13C12_MASSDIFF_U / precursorIon.charge;

        int smallestIsotope = std::max(0, int(std::ceil((lowerCutoff - precursorIon.monoMz) / isotopeStep)));
        int largestIsotope = int(std::floor((upperCutoff - precursorIon.monoMz) / isotopeStep));

        //loop through each isotope of precursor ion
        if (largestIsotope >= 0) {
            for (OpenMS::UInt m = smallestIsotope; m <= largestIsotope; ++m) {
                precursorIsotopes.insert(m);
            }
        }

        return precursorIsotopes;
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
    static void observedDistribution(std::vector<std::pair<double, double> > &obsDist,
                              std::vector<std::pair<double, double> > &theoDist,
                              const OpenMS::MSSpectrum<OpenMS::Peak1D> &spec)
    {
        //loop through each theoretical peak in isotopic distribution
        for (int i = 0; i < theoDist.size(); ++i) {

            //calculate search tolerance
            double tol =  OpenMS::Math::ppmToMass(ERROR_PPM, theoDist[i].first);

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
    std::vector<std::pair<double, double> > scaleDistribution(std::vector<std::pair<double, double> > &obsDist)
    {
        std::vector<std::pair<double, double> > scaled = obsDist;
        //sum intensities across distribution
        double totalIntensity = 0;
        for (int i = 0; i < scaled.size(); ++i) {
            totalIntensity += scaled[i].second;
        }

        //compute scaled intensity and replace value
        for (int j = 0; j < scaled.size(); ++j) {
            scaled[j].second = scaled[j].second / totalIntensity;
        }

        return scaled;
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
    static void exactConditionalFragmentIsotopeDist(std::vector<std::pair<double, double> > &condDist,
                                             const std::set<OpenMS::UInt> &precursorIsotopes,
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
            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / ion.charge ) * i;

            //set theoretical distribution pair
            std::pair<double, double> theo;
            theo.first = isoMZ;
            theo.second = condPeakList[i].second;
            condDist.push_back(theo);
        }
    }


    static void calibratedExactConditionalFragmentIsotopeDist(std::vector<std::pair<double, double> > &condDist,
                                                              OpenMS::Precursor &precursorInfo, double width,
                                                    const std::set<OpenMS::UInt> &precursorIsotopes,
                                                    const Ion &ion,
                                                    const Ion &precursorIon,
                                                              CalibrationModel &calModel)
    {
        //clear vector for distribution
        condDist.clear();


        OpenMS::UInt tMinIsotope = *std::min_element(precursorIsotopes.begin(), precursorIsotopes.end());
        OpenMS::UInt tMaxIsotope = *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end());

        OpenMS::UInt minIsotope = std::max(0, (int)tMinIsotope-3);
        OpenMS::UInt maxIsotope = tMaxIsotope+3;


        //compute conditional isotopic distribution and get vector of isotope peaks
        OpenMS::EmpiricalFormula precursorFormula = precursorIon.sequence.getFormula(OpenMS::Residue::Full, precursorIon.charge);

        double ionMZ = precursorIon.monoWeight / precursorIon.charge;

        for (int i = 0; i <= tMaxIsotope; ++i) {
            double fragIsoMZ = (ion.monoWeight / ion.charge) + ( OpenMS::Constants::C13C12_MASSDIFF_U / ion.charge ) * i;
            condDist.push_back(std::make_pair(fragIsoMZ, 0));
        }

        OpenMS::IsotopeDistribution prec_id = precursorIon.formula.getIsotopeDistribution(maxIsotope);
        prec_id.renormalize();

        for (int isotope = minIsotope; isotope <= maxIsotope; isotope++)
        {
            std::set<OpenMS::UInt> isotopes;
            isotopes.insert(isotope);
            OpenMS::IsotopeDistribution id = ion.formula.getConditionalFragmentIsotopeDist(precursorFormula, isotopes);
            id.renormalize();
            std::vector<std::pair<OpenMS::Size, double> > condPeakList = id.getContainer();

            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / precursorIon.charge ) * isotope;
            double offset = isoMZ - precursorInfo.getMZ();
            double efficiency = calModel.getEfficiency(width, offset) * prec_id.getContainer()[isotope].second;

            for (int i = 0; i < condPeakList.size() && i <= tMaxIsotope; ++i)
            {
                double calibratedIntensity = condPeakList[i].second * efficiency;
                condDist[i].second += calibratedIntensity;
            }
        }
    }

    static void approxPrecursorFromWeightIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                              const std::set<OpenMS::UInt> &precursorIsotopes,
                                              const Ion &fragmentIon)
    {
        OpenMS::UInt minIsotope = 7;
        //clear vector for distribution
        approxDist.clear();

        //construct distribution of depth at the maximum precursor isotope isolated
        OpenMS::IsotopeDistribution fragmentDist(std::min(minIsotope, *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) + 1));

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
            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / fragmentIon.charge ) * i;

            //set distribution pair
            std::pair<double, double> theo;
            theo.first = isoMZ;
            theo.second = isotopePeaks[i].second;
            approxDist.push_back(theo);
        }
    }

    static void approxFragmentFromWeightIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                             const std::set<OpenMS::UInt> &precursorIsotopes,
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
        OpenMS::IsotopeDistribution fragmentDist(*std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) + 1);

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
            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / fragmentIon.charge ) * i;

            //set distribution pair
            std::pair<double, double> theo;
            theo.first = isoMZ;
            theo.second = isotopePeaks[i].second;
            approxDist.push_back(theo);
        }
    }

    static void approxFragmentFromWeightAndSIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                                 const std::set<OpenMS::UInt> &precursorIsotopes,
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
        OpenMS::IsotopeDistribution fragmentDist(*std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) + 1);

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
            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / fragmentIon.charge ) * i;

            //set distribution pair
            std::pair<double, double> theo;
            theo.first = isoMZ;
            theo.second = isotopePeaks[i].second;
            approxDist.push_back(theo);
        }
    }

    static void approxFragmentSplineFromWeightIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                                    const std::set<OpenMS::UInt> &precursorIsotopes,
                                                    const Ion &fragmentIon,
                                                    const OpenMS::AASequence &precursorSequence,
                                                    const OpenMS::Int &precursorCharge,
                                                    const OpenMS::IsotopeSplineDB* isotopeDB)
    {
        //clear vector for distribution
        approxDist.clear();

        //precursor average weight
        double precursorAvgWeight = precursorSequence.getAverageWeight(OpenMS::Residue::Full, precursorCharge);
        //fragment average weight
        double fragmentAvgWeight = fragmentIon.formula.getAverageWeight();

        //construct distribution
        OpenMS::IsotopeDistribution fragmentDist(*std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) + 1);

        //estimate approx distribution from peptide weight
        fragmentDist = isotopeDB->estimateForFragmentFromPeptideWeight(precursorAvgWeight, fragmentAvgWeight, precursorIsotopes);
        //re-normalize distribution
        fragmentDist.renormalize();

        //get isotope vector
        std::vector<std::pair<OpenMS::Size, double> > isotopePeaks = fragmentDist.getContainer();

        //ion mz
        double ionMZ = fragmentIon.monoWeight / fragmentIon.charge;

        //loop through calculated isotopic distribution, fill with actual mz values
        for (int i = 0; i < isotopePeaks.size(); ++i) {

            //compute mz of isotope peak
            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / fragmentIon.charge ) * i;

            //set distribution pair
            std::pair<double, double> theo;
            theo.first = isoMZ;
            theo.second = isotopePeaks[i].second;
            approxDist.push_back(theo);
        }
    }

    static void approxFragmentSplineFromWeightAndSIsotopeDist(std::vector<std::pair<double, double> > &approxDist,
                                                        const std::set<OpenMS::UInt> &precursorIsotopes,
                                                        const Ion &fragmentIon,
                                                        const OpenMS::AASequence &precursorSequence,
                                                        const OpenMS::Int &precursorCharge,
                                                        const OpenMS::IsotopeSplineDB* isotopeDB)
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
        OpenMS::IsotopeDistribution fragmentDist(*std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) + 1);

        //estimate approx distribution from peptide weight
        fragmentDist = isotopeDB->estimateForFragmentFromPeptideWeightAndS(precursorAvgWeight, precursorSulfurs,
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
            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / fragmentIon.charge ) * i;

            //set distribution pair
            std::pair<double, double> theo;
            theo.first = isoMZ;
            theo.second = isotopePeaks[i].second;
            approxDist.push_back(theo);
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
    static void exactPrecursorIsotopeDist(std::vector<std::pair<double, double> > &theoDist,
                                          const std::set<OpenMS::UInt> &precursorIsotopes, const Ion &ion)
    {
        OpenMS::UInt minIsotope = 7;

        //construct distribution of depth at the maximum precursor isotope isolated
        OpenMS::UInt searchDepth = std::min(minIsotope, *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) + 1);

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
            double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / ion.charge ) * i;

            //set theoretical distribution pair
            std::pair<double, double> theo;
            theo.first = isoMZ;
            theo.second = theoPeakList[i].second;
            theoDist.push_back(theo);
        }
    }

    /**
     * Check if a scaled isotope distribution is following a characteristic isotope distribution.
     * @param dist scaled distribution where peak intensities sum to 1
     * @return true if distribution follows a characterist rise/fall
     */
    static bool scaledDistributionValid(const std::vector<std::pair<double, double> > &dist)
    {
        //return true;
        //distribution values decreasing flag
        bool valuesDecreasing = false;

        //loop from beginning of distribution, checking for decreasing values
        for (int i = 0; i < dist.size() - 1; ++i) {
            //if values are increasing
            if (!valuesDecreasing) {
                //check for next peak decreasing and difference is greater than 5%
                if ((dist[i+1].second < dist[i].second) && ((dist[i].second - dist[i+1].second) > 0.05) ) {
                    //next peak is greater than 5% less than current peak, distribution is decreasing
                    valuesDecreasing = true;
                }
            } else {    //values are decreasing
                //check for next peak increasing and difference is greater than 5%
                if ((dist[i+1].second > dist[i].second) && ((dist[i+1].second - dist[i].second) > 0.05) ) {
                    //distribution falling but next peak increases by more than 5%
                    return false;
                }
            }
        }
        //distribution follows normal rise and fall.
        return true;
    }
}


#endif //FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_SPECTRUMUTILITIES_H
