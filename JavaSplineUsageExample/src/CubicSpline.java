import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class CubicSpline {

    private List<Double> a_; // constant spline coefficients
    private List<Double> b_; // linear spline coefficients
    private List<Double> c_; // quadratic spline coefficients
    private List<Double> d_; // cubic spline coefficients
    private List<Double> x_; // knots

    public CubicSpline() {
        a_ = new ArrayList<>();
        b_ = new ArrayList<>();
        c_ = new ArrayList<>();
        d_ = new ArrayList<>();
        x_ = new ArrayList<>();
    }

    public CubicSpline(List<Double> a, List<Double> b, List<Double> c, List<Double> d, List<Double> x) {
        a_ = a;
        b_ = b;
        c_ = c;
        d_ = d;
        x_ = x;
    }

    public double eval(double x) {

        int index = Collections.binarySearch(x_, x);

        if (index < 0) {
            index = -(index+2);
        }

        double xx = x - x_.get(index);

        return ((d_.get(index) * xx + c_.get(index)) * xx + b_.get(index)) * xx + a_.get(index);
    }

    public boolean inBounds(double x) {
        return x >= x_.get(0) && x <= x_.get(x_.size()-1);
    }

}

