#include <fstream>
#include <iostream>
#include <sstream>
#include <OpenMS/FORMAT/MzMLFile.h>
#include "SpectrumUtilities.h"

struct TargetDetails {
    double center;
    double target;
    double offset;
    double width;

    TargetDetails() {};

    TargetDetails(double center, double target, double offset, double width)
            : center(center), target(target), offset(offset), width(width) {};
};

void usage() {

}

std::map<int, TargetDetails> parseInclusionList(std::string inclusionListPath) {
    std::map<int, TargetDetails> index2targetDetails;

    std::ifstream infile(inclusionListPath);

    std::string header;
    std::getline(infile, header);

    double center, target, offset, width;
    int z;
    std::string name, scanRange;

    int i = 0;
    while (infile >> center >> z >> name >> width >> scanRange) {

        int pos = name.find("_");
        target = std::stod(name.substr(0, pos));
        offset = std::stod(name.substr(pos+1, name.length()-pos));

        TargetDetails t(center, target, offset, width);
        index2targetDetails[i] = t;
        i++;
    }

    return index2targetDetails;
}

int main(int argc, char * argv[])
{

    std::map<int, TargetDetails> index2targetDetails = parseInclusionList(argv[2]);

    OpenMS::MzMLFile mzMLDataFile;
    OpenMS::MSExperiment<OpenMS::Peak1D> msExperiment;

    mzMLDataFile.load(argv[1], msExperiment);
    std::ofstream out(argv[3]);

    double ppmTol = std::atof(argv[4]);

    out << "intensity" << "\t" << "target" << "\t" << "offset" << "\t" << "width" << std::endl;

    for (int specIndex = 0; specIndex < msExperiment.getNrSpectra(); ++specIndex) {
        //get copy of current spectrum
        OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrum = msExperiment.getSpectrum(specIndex);

        //sort spectrum by mz
        currentSpectrum.sortByPosition();

        TargetDetails targetDetails = index2targetDetails[specIndex % index2targetDetails.size()];

        double tol = OpenMS::Math::ppmToMass(ppmTol, targetDetails.target);

        //find nearest peak to ion mz within tolerance
        OpenMS::Int peakIndex = currentSpectrum.findNearest(targetDetails.target, tol);

        double intensity = peakIndex == -1 ? 0 : currentSpectrum[peakIndex].getIntensity();
        out << intensity << "\t" << targetDetails.target << "\t" << targetDetails.offset << "\t" << targetDetails.width << std::endl;
    }

    std::cout << "done" << std::endl;
}