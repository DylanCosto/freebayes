#ifndef FREEBAYES_EWENS_H
#define FREEBAYES_EWENS_H

#include <map>
#include <cmath>
#include "Utility.h"

using namespace std;

// genotype priors

long double alleleFrequencyProbability(const map<int, int>& alleleFrequencyCounts, long double theta);
long double alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta);
long double impl_alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta);

// Keyed by theta as well as the counts, so that a cached value is never used
// for a different theta
class AlleleFrequencyProbabilityCache : public map<long double, map<map<int, int>, long double> > {
public:
    long double alleleFrequencyProbabilityln(const map<int, int>& counts, long double theta) {
        map<map<int, int>, long double>& t = (*this)[theta];
        map<map<int, int>, long double>::iterator p = t.find(counts);
        if (p == t.end()) {
            long double pln = impl_alleleFrequencyProbabilityln(counts, theta);
            t.insert(make_pair(counts, pln));
            return pln;
        } else {
            return p->second;
        }
    }
};

#endif
