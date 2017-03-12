#include <iostream>
#include <chrono>
#include <random>
#include <vector>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis;

static const OpenMS::IsotopeSplineDB* splineDB = OpenMS::IsotopeSplineDB::getInstance();

void timePrecursorSpline(std::vector<double> &masses, UInt max_depth, bool print)
{
    for (UInt depth = 1; depth <= max_depth; ++depth)
    {
        auto time_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < masses.size()-1; ++i) {
            IsotopeDistribution id = splineDB->estimateFromPeptideWeight(masses[i], depth);
        }
        auto time_end = std::chrono::high_resolution_clock::now();

        std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start);

        if (print) {
            std::cout << "Spline" << "\t" << d.count() << "\t" << depth << "\t" << "Precursor masses" << std::endl;
        }
    }
}

void timePrecursorFFT(std::vector<double> &masses, UInt max_depth, bool print)
{
    for (UInt depth = 1; depth <= max_depth; ++depth)
    {
        auto time_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < masses.size()-1; ++i) {
            IsotopeDistribution id(depth);
            id.estimateFromPeptideWeight(masses[i]);
        }
        auto time_end = std::chrono::high_resolution_clock::now();

        std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start);

        if (print) {
            std::cout << "FFT" << "\t" << d.count() << "\t" << depth << "\t" << "Precursor masses" << std::endl;
        }
    }
}

void timeFragmentFFTSingle(std::vector<double> &masses, UInt max_depth, bool print)
{
    std::set<UInt> precursor_isotopes;
    for (UInt i = 0; i < max_depth; ++i) {
        precursor_isotopes.clear();
        precursor_isotopes.insert(i);
        auto time_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < masses.size()-1; ++i) {
            IsotopeDistribution id;
            id.estimateForFragmentFromPeptideWeight(masses[i]+masses[i+1], masses[i], precursor_isotopes);
        }
        auto time_end = std::chrono::high_resolution_clock::now();

        std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start);

        if (print) {
            std::cout << "FFT" << "\t" << d.count() << "\t" << i+1 << "\t" << "Single fragment isotope" << std::endl;
        }
    }
}

void timeFragmentFFTCombined(std::vector<double> &masses, UInt max_depth, bool print)
{
    std::set<UInt> precursor_isotopes;
    for (UInt i = 0; i < max_depth; ++i)
    {
        precursor_isotopes.insert(i);
        auto time_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < masses.size()-1; ++i) {
            IsotopeDistribution id;
            id.estimateForFragmentFromPeptideWeight(masses[i]+masses[i+1], masses[i], precursor_isotopes);
        }
        auto time_end = std::chrono::high_resolution_clock::now();

        std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start);

        if (print) {
            std::cout << "FFT" << "\t" << d.count() << "\t" << i+1 << "\t" << "Multiple fragment isotopes" << std::endl;
        }
    }
}

void timeFragmentSplineSingle(std::vector<double> &masses, UInt max_depth, bool print) {
    std::set<UInt> precursor_isotopes;
    for (UInt i = 0; i < max_depth; ++i) {
        precursor_isotopes.clear();
        precursor_isotopes.insert(i);
        auto time_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < masses.size()-1; ++i) {
            IsotopeDistribution id = splineDB->estimateForFragmentFromPeptideWeight(masses[i]+masses[i+1], masses[i], precursor_isotopes);
        }
        auto time_end = std::chrono::high_resolution_clock::now();

        std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start);

        if (print) {
            std::cout << "Spline" << "\t" << d.count() << "\t" << i+1 << "\t" << "Single fragment isotope" << std::endl;
        }
    }
}

void timeFragmentSplineCombined(std::vector<double> &masses, UInt max_depth, bool print)
{
    std::set<UInt> precursor_isotopes;
    for (UInt i = 0; i < max_depth; ++i)
    {
        precursor_isotopes.insert(i);
        auto time_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < masses.size()-1; ++i) {
            IsotopeDistribution id = splineDB->estimateForFragmentFromPeptideWeight(masses[i]+masses[i+1], masses[i], precursor_isotopes);
        }
        auto time_end = std::chrono::high_resolution_clock::now();

        std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start);

        if (print) {
            std::cout << "Spline" << "\t" << d.count() << "\t" << i+1 << "\t" << "Multiple fragment isotopes" << std::endl;
        }
    }


}

void usage()
{

}

std::vector<double> getRandomMasses(int num_tests) {
    std::vector<double> masses;
    for (int i = 0; i < num_tests; ++i) {
        masses.push_back(dis(gen));
    }
    return masses;
}

int main(int argc, const char ** argv) {

    float min_mass = atof(argv[1]);
    float max_mass = atof(argv[2]);
    UInt max_depth = atoi(argv[3]);
    int num_tests =  atoi(argv[4]);

    dis = std::uniform_real_distribution<>(min_mass, max_mass);

    // Print header
    std::cout << "method\t" << "time\t" << "isotope\t" << "comparison" << std::endl;

    std::vector<double> masses = getRandomMasses(num_tests);
    // Time fragment FFT Single
    timeFragmentFFTCombined(masses, max_depth, false);
    timeFragmentFFTCombined(masses, max_depth, true);
    // Time fragment spline Combined
    timeFragmentSplineCombined(masses, max_depth, false);
    timeFragmentSplineCombined(masses, max_depth, true);


    // Time fragment FFT Combined
    timeFragmentFFTSingle(masses, max_depth, false);
    timeFragmentFFTSingle(masses, max_depth, true);
    // Time fragment spline Single
    timeFragmentSplineSingle(masses, max_depth, false);
    timeFragmentSplineSingle(masses, max_depth, true);


    // Time precursor FFT
    timePrecursorFFT(masses, max_depth, false);
    timePrecursorFFT(masses, max_depth, true);
    // Time precursor spline
    timePrecursorSpline(masses, max_depth, false);
    timePrecursorSpline(masses, max_depth, true);

    return 0;
}