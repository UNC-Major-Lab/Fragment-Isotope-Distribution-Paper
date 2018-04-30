import java.util.List;
import java.util.Map;

public class IsotopeSplineDB {

    private static final double C13C12_MASSDIFF_U = 1.0033548378;

    // a key of -1 is for the average splines
    // the index in List<CublicSpline> is the isotope
    private Map<Integer, List<CubicSpline>> numSulfur2models_;

    public IsotopeSplineDB(String splinePath) {
        readSplinesFromFile_(splinePath);
    }

    private void readSplinesFromFile_(String splinePath) {
        IsotopeSplineXMLParser parser = new IsotopeSplineXMLParser();

        numSulfur2models_ = parser.parse(splinePath);
    }

    IsotopeDistribution estimateFromPeptideWeight(double monoMass, int maxIsotope) {
        return estimateFromPeptideWeightAndSulfur(monoMass, maxIsotope, -1);
    }

    IsotopeDistribution estimateFromPeptideWeightAndSulfur(double monoMass, int maxIsotope, int numSulfur) {
        IsotopeDistribution result;

        if (inModelBounds(monoMass, maxIsotope, numSulfur)) {
            result = new IsotopeDistribution(maxIsotope);
            for (int isotope = 0; isotope <= maxIsotope; ++isotope) {
                double probability = numSulfur2models_.get(numSulfur).get(isotope).eval(monoMass);
                result.setMass(isotope, monoMass + (isotope * C13C12_MASSDIFF_U));
                result.setIntensity(isotope, (float) probability);
            }
        } else {
            throw new RuntimeException("Request out of range.");
        }

        return result;
    }

    public IsotopeDistribution estimateForFragmentFromWeights(double monoPeptideMass, double monoFragmentMass,
                                                              int minIsotope, int maxIsotope) {

        return estimateForFragmentFromWeightsAndSulfur(monoPeptideMass, monoFragmentMass, minIsotope, maxIsotope, -1, -1);
    }

    public IsotopeDistribution estimateForFragmentFromWeightsAndSulfur(double monoPeptideMass, double monoFragmentMass,
                                                              int minIsotope, int maxIsotope, int precursorSulfur, int fragmentSulfur) {

        IsotopeDistribution fragment = estimateFromPeptideWeightAndSulfur(monoFragmentMass, maxIsotope, precursorSulfur);
        IsotopeDistribution compFragment = estimateFromPeptideWeightAndSulfur(monoPeptideMass - monoFragmentMass, maxIsotope, fragmentSulfur);
        return calcFragmentIsotopeDistribution(fragment, compFragment, minIsotope, maxIsotope);

    }

    public boolean inModelBounds(double monoMass, int maxIsotope, int numSulfur) {
        // Check if we have a sulfur-specific model for this
        if (!numSulfur2models_.containsKey(numSulfur)) {
            return false;
        }

        // Check if max isotope is in bounds
        if (maxIsotope > numSulfur2models_.get(numSulfur).size()-1) {
            return false;
        }

        // Check if masses are in bounds
        for (int isotope = 0; isotope <= maxIsotope; ++isotope) {
            if (!numSulfur2models_.get(numSulfur).get(isotope).inBounds(monoMass)) {
                return false;
            }
        }

        // All checks passed
        return true;
    }

    private IsotopeDistribution calcFragmentIsotopeDistribution(IsotopeDistribution fragment, IsotopeDistribution compFragment,
                                                                int minIsotope, int maxIsotope) {

        IsotopeDistribution result = new IsotopeDistribution(maxIsotope);

        for (int i = 0; i < fragment.size(); ++i) {
            for (int isotope = minIsotope; isotope <= maxIsotope; ++isotope) {
                if (isotope >= i && (isotope-i) < compFragment.size()) {
                    result.setIntensity(i, (float) (result.getIntensity(i) + compFragment.getIntensity(isotope-i)));
                }
            }
            result.setIntensity(i, (float) (result.getIntensity(i) * fragment.getIntensity(i)));
            result.setMass(i, fragment.getMass(0) + (i * C13C12_MASSDIFF_U));
        }

        return result;
    }

}