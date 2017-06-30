//
// Created by Dennis Goldfarb on 6/9/17.
//

#ifndef FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_ISOTOPEDISTRIBUTIONS_H
#define FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_ISOTOPEDISTRIBUTIONS_H


#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include "SpectrumUtilities.h"
#include "Stats.h"
#include "CalibrationModel.h"

class IsotopeDistributions {

public :

    // Precursor distributions
    IsotopeDistributions(std::set<OpenMS::UInt> precursorIsotopes, Ion precursorIon,
                         const OpenMS::IsotopeSplineDB* isotopeDB,
                         OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumCentroid,
                         CalibrationModel &calModel, double width, OpenMS::Precursor &precursorInfo)
    {
        if (precursorIsotopes.size() > 0) {
            OpenMS::UInt minIsotope = *std::min_element(precursorIsotopes.begin(), precursorIsotopes.end());
            OpenMS::UInt searchDepth = *std::max_element(precursorIsotopes.begin(), precursorIsotopes.end()) + 1;

            OpenMS::IsotopeDistribution id = precursorIon.formula.getIsotopeDistribution(searchDepth);
            double ionMZ = precursorIon.monoWeight / precursorIon.charge;

            for (int i = minIsotope; i < id.size(); ++i) {
                double isoMZ = ionMZ + ( OpenMS::Constants::C13C12_MASSDIFF_U / precursorIon.charge ) * i;
                exactPrecursorDist.push_back(std::make_pair(isoMZ, id.getContainer()[i].second));

                double offset = isoMZ - precursorInfo.getMZ();
                double efficiency = calModel.getEfficiency(width, offset);
                double calibratedIntensity = exactPrecursorDist.back().second * efficiency;
                calibratedExactPrecursorDist.push_back(std::make_pair(isoMZ, calibratedIntensity));
            }
            exactPrecursorDist = SpectrumUtilities::scaleDistribution(exactPrecursorDist);
            calibratedExactPrecursorDist = SpectrumUtilities::scaleDistribution(calibratedExactPrecursorDist);

            SpectrumUtilities::observedDistribution(observedDist, exactPrecursorDist, currentSpectrumCentroid);
            scaledObservedDist = SpectrumUtilities::scaleDistribution(observedDist);
        }

        exactPrecursorX2 = Stats::computeX2(scaledObservedDist, exactPrecursorDist);
        calibratedExactPrecursorX2 = Stats::computeX2(scaledObservedDist, calibratedExactPrecursorDist);
    }


    // Fragment distributions
    IsotopeDistributions(std::set<OpenMS::UInt> precursorIsotopes, Ion ion, Ion precursorIon,
                         const OpenMS::IsotopeSplineDB* isotopeDB, OpenMS::MSSpectrum<OpenMS::Peak1D> &currentSpectrumCentroid,
                         CalibrationModel &calModel, OpenMS::Precursor &precursorInfo, double width)
    {

        SpectrumUtilities::exactConditionalFragmentIsotopeDist(exactConditionalFragmentDist,
                                                               precursorIsotopes,
                                                               ion,
                                                               precursorIon.sequence,
                                                               precursorIon.charge);

        SpectrumUtilities::approxFragmentFromWeightIsotopeDist(approxFragmentFromWeightDist,
                                                               precursorIsotopes,
                                                               ion,
                                                               precursorIon.sequence,
                                                               precursorIon.charge);

        SpectrumUtilities::approxFragmentFromWeightAndSIsotopeDist(approxFragmentFromWeightAndSulfurDist,
                                                                   precursorIsotopes,
                                                                   ion,
                                                                   precursorIon.sequence,
                                                                   precursorIon.charge);

        SpectrumUtilities::approxFragmentSplineFromWeightIsotopeDist(approxFragmentSplineFromWeightDist,
                                                                     precursorIsotopes,
                                                                     ion,
                                                                     precursorIon.sequence,
                                                                     precursorIon.charge,
                                                                     isotopeDB);

        SpectrumUtilities::approxFragmentSplineFromWeightAndSIsotopeDist(approxFragmentSplineFromWeightAndSulfurDist,
                                                                         precursorIsotopes,
                                                                         ion,
                                                                         precursorIon.sequence,
                                                                         precursorIon.charge,
                                                                         isotopeDB);

        SpectrumUtilities::exactPrecursorIsotopeDist(exactPrecursorDist, precursorIsotopes,
                                                     ion);

        SpectrumUtilities::approxPrecursorFromWeightIsotopeDist(approxPrecursorFromWeightDist, precursorIsotopes,
                                                                ion);

        //match theoretical distribution with observed peaks
        SpectrumUtilities::observedDistribution(observedDist, exactConditionalFragmentDist, currentSpectrumCentroid);
        //scale observed intensities across distribution
        scaledObservedDist = SpectrumUtilities::scaleDistribution(observedDist);


        //report on matched ion distribution depth
        completeFlag = true;
        completeAtDepth = 0;
        for (int i = 0; i < observedDist.size(); ++i) {
            if (observedDist[i].second != 0) {
                if (completeFlag) {
                    ++completeAtDepth;
                }
            } else {
                completeFlag = false;
            }
        }

        isValid = SpectrumUtilities::scaledDistributionValid(scaledObservedDist) && completeAtDepth > 1;

        //compute chi-squared with observed to Conditional
        exactCondFragmentX2 = Stats::computeX2(scaledObservedDist, exactConditionalFragmentDist);
        approxFragmentFromWeightX2 = Stats::computeX2(scaledObservedDist, approxFragmentFromWeightDist);
        approxFragmentFromWeightAndSulfurX2 = Stats::computeX2(scaledObservedDist, approxFragmentFromWeightAndSulfurDist);
        approxFragmentSplineFromWeightX2 = Stats::computeX2(scaledObservedDist, approxFragmentSplineFromWeightDist);
        approxFragmentSplineFromWeightAndSulfurX2 = Stats::computeX2(scaledObservedDist, approxFragmentSplineFromWeightAndSulfurDist);
        exactPrecursorX2 = Stats::computeX2(scaledObservedDist, exactPrecursorDist);
        approxPrecursorX2 = Stats::computeX2(scaledObservedDist, approxPrecursorFromWeightDist);

        if (calModel.width2offset2efficiency.size() > 0) {
            SpectrumUtilities::calibratedExactConditionalFragmentIsotopeDist(calibratedExactConditionalFragmentDist,
                                                                             precursorInfo, width,
                                                                             precursorIsotopes,
                                                                             ion,
                                                                             precursorIon,
                                                                             calModel);
            calibratedExactConditionalFragmentDist = SpectrumUtilities::scaleDistribution(
                    calibratedExactConditionalFragmentDist);
            calibratedExactFragmentX2 = Stats::computeX2(scaledObservedDist, calibratedExactConditionalFragmentDist);
        }



    }

    // Precursor or fragment
    std::vector<std::pair<double, double> > exactPrecursorDist;
    std::vector<std::pair<double, double> > approxPrecursorFromWeightDist;
    std::vector<std::pair<double, double> > observedDist;
    std::vector<std::pair<double, double> > scaledObservedDist;

    // Precursor
    std::vector<std::pair<double, double> > approxPrecursorFromWeightAndSulfurDist;
    std::vector<std::pair<double, double> > approxPrecursorSplineFromWeightDist;
    std::vector<std::pair<double, double> > approxPrecursorSplineFromWeightAndSulfurDist;
    std::vector<std::pair<double, double> > calibratedExactPrecursorDist;

    // Fragment
    std::vector<std::pair<double, double> > approxFragmentFromWeightDist;
    std::vector<std::pair<double, double> > approxFragmentFromWeightAndSulfurDist;
    std::vector<std::pair<double, double> > approxFragmentSplineFromWeightDist;
    std::vector<std::pair<double, double> > approxFragmentSplineFromWeightAndSulfurDist;
    std::vector<std::pair<double, double> > exactConditionalFragmentDist;
    std::vector<std::pair<double, double> > calibratedExactConditionalFragmentDist;


    // Precursor or fragment
    double exactPrecursorX2;
    double approxPrecursorX2;

    // Precursor
    double approxPrecursorAndSulfurX2;
    double approxPrecursorSplineX2;
    double approxPrecursorSplineAndSulfurX2;
    double calibratedExactPrecursorX2;

    // Fragment
    double exactCondFragmentX2;
    double approxFragmentFromWeightX2;
    double approxFragmentFromWeightAndSulfurX2;
    double approxFragmentSplineFromWeightX2;
    double approxFragmentSplineFromWeightAndSulfurX2;
    double calibratedExactFragmentX2;


    bool completeFlag;
    int completeAtDepth;

    bool isValid;


private:

};


#endif //FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_ISOTOPEDISTRIBUTIONS_H
