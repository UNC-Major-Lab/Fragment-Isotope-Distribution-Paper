public class Main {

    //public static final IsotopeSplineDB splineDB = new IsotopeSplineDB("data/IsotopeSplines_10kDa_21isotopes.xml");
    public static final IsotopeSplineDB splineDB = new IsotopeSplineDB("data/IsotopeSplines_100kDa_101isotopes.xml");

    public static void main(String [] args) {

        IsotopeDistribution id;

        id = splineDB.estimateFromPeptideWeight(1000, 5);
        System.out.println("Approximate precursor distribution for mass: " + 1000 + " and maxIsotope: " + 5);
        printDistribution(id);

        id = splineDB.estimateFromPeptideWeightAndSulfur(1000, 5, 0 );
        System.out.println("Approximate precursor distribution for mass: " + 1000 + " and maxIsotope: " + 5 + " and numSulfurs: " + 0);
        printDistribution(id);

        id = splineDB.estimateForFragmentFromWeights(1000, 500, 1, 2 );
        System.out.println("Approximate fragment distribution for precursorMass: " + 1000 + " and fragmentMass: " + 5 +
                " and minPrecursorIsotope: " + 1 + " and maxPrecursorIsotope: " + 2);
        printDistribution(id);

        id = splineDB.estimateForFragmentFromWeightsAndSulfur(1000, 500, 1, 2, 2, 0 );
        System.out.println("Approximate fragment distribution for precursorMass: " + 1000 + " and fragmentMass: " + 5 +
                " and minPrecursorIsotope: " + 1 + " and maxPrecursorIsotope: " + 2 + " and precursorNumSulfurs: " + 2 +
                " and fragmentNumSulfurs: " + 0);
        printDistribution(id);

        // 50kDa will only work with the spline file trained on data up to 100kDa
        id = splineDB.estimateFromPeptideWeight(50000, 50);
        System.out.println("Approximate precursor distribution for mass: " + 50000 + " and maxIsotope: " + 50);
        printDistribution(id);
    }


    static public void printDistribution(IsotopeDistribution id) {
        for (int i = 0; i < id.size(); ++i) {
            System.out.println(id.getMass(i) + "\t" + id.getIntensity(i));
        }
    }


}
