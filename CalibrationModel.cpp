//
// Created by Dennis Goldfarb on 6/14/17.
//

#include "CalibrationModel.h"

#include <cmath>
#include <numeric>

void CalibrationModel::finishCalibation()
{
    for (auto &width_itr : width2offset2intensities)
    {
        double maxIntensity = 0.0;

        for (auto &offset_itr : width_itr.second)
        {
            std::vector<double> intensities = offset_itr.second;

            std::sort(intensities.begin(), intensities.end());

            std::vector<double> keepers = intensities;

            //remove outliers
            /*int N = intensities.size();
            double mean = std::accumulate(intensities.begin(), intensities.end(), 0.0)/ N;
            double sq_sum = std::inner_product(intensities.begin(), intensities.end(), intensities.begin(), 0.0);
            double stdev = std::sqrt(sq_sum / intensities.size() - mean * mean);


            for (int i = 0; i < N; ++i) {
                double z = (intensities[i] - mean) / stdev;
                if (z <= 3.5) {
                    keepers.push_back(intensities[i]);
                } else {
                    int x = 1;
                }
            }*/


            //median
            if (keepers.size() > 0) {
                if (keepers.size() % 2 == 1) {
                    width2offset2efficiency[width_itr.first][offset_itr.first] = keepers[keepers.size() / 2];
                } else {
                    width2offset2efficiency[width_itr.first][offset_itr.first] =
                            (keepers[keepers.size() / 2] + keepers[1 + keepers.size() / 2]) / 2;
                }

                if (width2offset2efficiency[width_itr.first][offset_itr.first] > maxIntensity) {
                    maxIntensity = width2offset2efficiency[width_itr.first][offset_itr.first];
                }
            }
        }

        for (auto &offset_itr : width_itr.second)
        {
            width2offset2efficiency[width_itr.first][offset_itr.first] /= maxIntensity;
        }
    }


}

void CalibrationModel::addCalibrationPoint(double width, double offset, double intensity)
{
    width = std::round(width * 100) / 100;
    width2offset2intensities[width][offset].push_back(intensity);
}

double CalibrationModel::getEfficiency(double width, double offset) {
    width = std::round(width * 100) / 100;

    // find nearest offsets and linearly interpolate in between
    double leftOffset, rightOffset, leftDiff = -1000, rightDiff = 1000;
    for (auto offset_itr : width2offset2efficiency[width])
    {
        double newDiff = offset_itr.first - offset;

        if (newDiff < 0 && newDiff > leftDiff)
        {
            leftDiff = newDiff;
            leftOffset = offset_itr.first;
        }

        if (newDiff > 0 && newDiff < rightDiff)
        {
            rightDiff = newDiff;
            rightOffset = offset_itr.first;
        }
    }

    if (leftDiff == -1000 || rightDiff == 1000) return 0;

    double totalDiff = rightOffset - leftOffset;
    double percentLeftDiff = 1 - (std::abs(leftDiff) / totalDiff);
    double percentRightDiff = 1 - (std::abs(rightDiff) / totalDiff);

    double efficiency = (width2offset2efficiency[width][leftOffset] * percentLeftDiff)
                         + (width2offset2efficiency[width][rightOffset] * percentRightDiff);

    return efficiency;
}
