//
// Created by Dennis Goldfarb on 9/16/16.
//

#ifndef EXAMPLE_PROJECT_USING_OPENMS_STATS_H
#define EXAMPLE_PROJECT_USING_OPENMS_STATS_H

#include <cmath>

class Stats {

public:
    template <typename IteratorType1, typename IteratorType2>
    static double totalVariationDistance(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
    {
        double sum = 0;
        IteratorType1 iter_a = begin_a;
        IteratorType2 iter_b = begin_b;
        for (; iter_a != end_a; ++iter_a, ++iter_b)
        {
            sum += std::abs(*iter_a - *iter_b);
        }
        return sum;
    }

    template <typename IteratorType1, typename IteratorType2>
    static double chiSquared(IteratorType1 begin_o, IteratorType1 end_o, IteratorType2 begin_e, IteratorType2 end_e)
    {
        double sum = 0;
        IteratorType1 iter_o = begin_o;
        IteratorType2 iter_e = begin_e;
        for (; iter_o != end_o; ++iter_o, ++iter_e)
        {
            if (*iter_e != 0)
            {
                double diff = *iter_o - *iter_e;
                sum += (diff * diff) / *iter_e;
            }
        }
        return sum;
    }

    static double computeCC(const std::vector<std::pair<double, double> > &obsDist,
                     const std::vector<std::pair<double, double> > &theoDist)
    {
        //vector to hold observed proportions
        std::vector<double> obsProp;
        //vector to hold theoretical proportions
        std::vector<double> theoProp;

        //check they are both the same size
        if (obsDist.size() != theoDist.size()) {
            return 0;
        }

        //fill proportions vectors from distribution parameters
        for (int i = 0; i < obsDist.size(); ++i) {
            obsProp.push_back(obsDist[i].second);
            theoProp.push_back(theoDist[i].second);
        }

        //compute pearsons correlation coefficient
        return OpenMS::Math::pearsonCorrelationCoefficient(obsProp.begin(), obsProp.end(),
                                                           theoProp.begin(), theoProp.end());
    }

    static double computeX2(const std::vector<std::pair<double, double> > &obsDist,
                     const std::vector<std::pair<double, double> > &theoDist)
    {
        //vector to hold observed proportions
        std::vector<double> obsProp;
        //vector to hold theoretical proportions
        std::vector<double> theoProp;

        //check they are both the same size
        if (obsDist.size() != theoDist.size()) {
            return -1;
        }

        //fill proportions vectors from distribution parameters
        for (int i = 0; i < obsDist.size(); ++i) {
            obsProp.push_back(obsDist[i].second);
            theoProp.push_back(theoDist[i].second);
        }

        //compute chi squared statistic
        return Stats::chiSquared(obsProp.begin(), obsProp.end(),
                                 theoProp.begin(), theoProp.end());
    }

    static double computeVD(const std::vector<std::pair<double, double> > &obsDist,
                     const std::vector<std::pair<double, double> > &theoDist)
    {
        //vector to hold observed proportions
        std::vector<double> obsProp;
        //vector to hold theoretical proportions
        std::vector<double> theoProp;

        //check they are both the same size
        if (obsDist.size() != theoDist.size()) {
            return -1;
        }

        //fill proportions vectors from distribution parameters
        for (int i = 0; i < obsDist.size(); ++i) {
            obsProp.push_back(obsDist[i].second);
            theoProp.push_back(theoDist[i].second);
        }

        //compute total variation distance
        return Stats::totalVariationDistance(obsProp.begin(), obsProp.end(),
                                             theoProp.begin(), theoProp.end());
    }
};


#endif //EXAMPLE_PROJECT_USING_OPENMS_STATS_H
