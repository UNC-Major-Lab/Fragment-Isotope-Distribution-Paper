#include <iostream>
#include <chrono>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

using namespace OpenMS;

void testFFTApproximationGivenIsotopes(float min_mass, float max_mass, std::vector<UInt> precursor_isotopes)
{
    for (float precursor_mass = min_mass + 1; precursor_mass < max_mass; precursor_mass+=100)
    {
        for (float fragment_mass = min_mass; fragment_mass < precursor_mass; fragment_mass+=100)
        {
            IsotopeDistribution id;
            id.estimateForFragmentFromPeptideWeightFast(precursor_mass, fragment_mass, precursor_isotopes);
        }
    }
}

void testFFTApproximation(float min_mass, float max_mass, UInt max_isotope)
{
    std::vector<UInt> precursor_isotopes;
    for (UInt i = 0; i <= max_isotope; ++i)
    {
        precursor_isotopes.clear();
        precursor_isotopes.push_back(i);
        testFFTApproximationGivenIsotopes(min_mass, max_mass, precursor_isotopes);
    }

    precursor_isotopes.clear();
    for (UInt i = 0; i <= max_isotope; ++i)
    {
        precursor_isotopes.push_back(i);
        testFFTApproximationGivenIsotopes(min_mass, max_mass, precursor_isotopes);
    }
}

void usage()
{

}

int main(int argc, const char ** argv) {

    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::milliseconds ms;


    auto time_start = Time::now();

    testFFTApproximation(200, 8500, 0);

    auto time_end = Time::now();

    std::chrono::duration<float> duration = time_end - time_start;
    ms d = std::chrono::duration_cast<ms>(duration);

    std::cout << duration.count() << "s" << std::endl;
    std::cout << d.count() << "ms" << std::endl;


    return 0;
}