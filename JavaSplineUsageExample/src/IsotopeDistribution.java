import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class IsotopeDistribution {

    private List<Double> massData;
    private List<Double> intensityData;

    public IsotopeDistribution() {
        massData = new ArrayList();
        intensityData = new ArrayList();
    }

    IsotopeDistribution(int maxIsotope) {
        massData = new ArrayList(maxIsotope+1);
        intensityData = new ArrayList(maxIsotope+1);
        for (int i = 0; i <= maxIsotope; ++i) {
            massData.add(0.0);
            intensityData.add(0.0);
        }
    }

    public int size() {
        return massData.size();
    }

    public double getIntensity(int i) {
        return intensityData.get(i);
    }

    void setIntensity(int i, double intensity) {
        intensityData.set(i, intensity);
    }

    void setMass(int i, double mz) {
        massData.set(i, mz);
    }

    public double getMass(int i) {
        return massData.get(i);
    }

    public List<Double> getMassData() {
        return massData;
    }

    public List<Double> getIntensityData() {
        return intensityData;
    }

    public void normalizeToBasePeak() {
        normalizeToValue(Collections.max(intensityData));
    }

    public void normalizeToValue(double value) {
        for (int i = 0; i < intensityData.size(); ++i) {
            intensityData.set(i, intensityData.get(i) / value);
        }
    }
}