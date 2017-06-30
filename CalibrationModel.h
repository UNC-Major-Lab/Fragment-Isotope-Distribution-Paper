//
// Created by Dennis Goldfarb on 6/14/17.
//

#ifndef FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_CALIBRATIONMODEL_H
#define FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_CALIBRATIONMODEL_H


#include <map>
#include <vector>

class CalibrationModel {

public:
    CalibrationModel() {}

    double getEfficiency(double width, double offset);

    void addCalibrationPoint(double width, double offset, double intensity);

    void finishCalibation();

    std::map<double, std::map<double, double> > width2offset2efficiency;
    std::map<double, std::map<double, std::vector<double> > > width2offset2intensities;
};


#endif //FRAGMENT_ISOTOPE_DISTRIBUTION_PAPER_CALIBRATIONMODEL_H
