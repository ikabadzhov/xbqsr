
#include "TablesBQSR.h"

const int32_t LOG10_QUAL_EMPIRICAL_PRIOR_CACHE_SIZE = 21;

const double LOG10_QUAL_EMPIRICAL_PRIOR_CACHE[LOG10_QUAL_EMPIRICAL_PRIOR_CACHE_SIZE] = {
            -0.045757490560675115,
            -0.9143464543671788,
            -3.5201133457866898,
            -7.863058164819208,
            -13.943180911464733,
            -21.760481585723266,
            -31.314960187594806,
            -42.606616717079355,
            -55.63545117417691,
            -70.40146355888747,
            -86.90465387121104,
            -105.14502211114761,
            -125.1225682786972,
            -146.83729237385978,
            -170.2891943966354,
            -195.47827434702398,
            -222.4045322250256,
            -251.06796803064023,
            -281.46858176386786,
            -313.60637342472336,
            -1.7976931348623157e308,
};

//Log10QualEmpiricalPriorCache log10QualEmpiricalPriorCache;

bool bqsrTableKey::operator==(const bqsrTableKey &rhs) const {
    return Qual == rhs.Qual &&
           Covariate == rhs.Covariate &&
           ReadGroup == rhs.ReadGroup;
}

bool bqsrTableKey::operator!=(const bqsrTableKey &rhs) const {
    return !(rhs == *this);
}


/**
 * bqsrTable functions
 * */

void bqsrTable::update(bqsrEntry* table, size_t idx, int64_t errorVal){
    if (table[idx].Observations == -1) {
        table[idx].Observations++;
    }
    table[idx].Observations++;
    table[idx].Mismatches += errorVal;
}

/*size_t bqsrTable::size() {
    return table.size();
}*/


uint8_t  calculateBayesianEstimateOfEmpiricalQuality(int64_t observations, int64_t mismatches, double priorQuality) {
    if (observations > maxNumberOfObservations) {
        mismatches = static_cast<int64_t>(round(static_cast<double>(mismatches) * (static_cast<double>(maxNumberOfObservations) / static_cast<double>(observations))));
        observations = maxNumberOfObservations;
    }
    volatile uint8_t nbins = maxReasonableQualScore + 1;
    volatile double max = -std::numeric_limits<double>::max();
    uint8_t maxI = 0;
    for (uint8_t i = 0; i < nbins; i++) {
        auto floati = static_cast<double>(i);
        double p1 = log10QualEmpiricalPrior(floati, priorQuality);
        double p2 = log10QualEmpiricalLikelihood(floati, observations, mismatches);
        double log10Posterior = p1 + p2;
        if (max < log10Posterior) {
            max = log10Posterior;
            maxI = i;
        }
    }
    return maxI;
}

uint8_t  bqsrEntry::calculateEmpiricalQuality(double priorQual) {
    int64_t smoothMismatches = this->Mismatches + smoothingConstant;
    int64_t smoothObservations = this->Observations + smoothingConstant + smoothingConstant;
    uint8_t empiricalQual = calculateBayesianEstimateOfEmpiricalQuality(smoothObservations, smoothMismatches, priorQual);
    return std::min(empiricalQual, static_cast<uint8_t>(maxRecalibratedQualScore));
}

double log10QualEmpiricalPrior(double empiricalQual, double reportedQual) {
    auto diff = std::min(std::abs(static_cast<int32_t>(empiricalQual-reportedQual)), static_cast<int32_t>(LOG10_QUAL_EMPIRICAL_PRIOR_CACHE_SIZE-1));
    return LOG10_QUAL_EMPIRICAL_PRIOR_CACHE[diff];
}

double log10Gamma(int64_t n) {
    double gammaN = lgamma(static_cast<double>(n));
    return gammaN * Log10E;
}

double log10BinomialCoefficient(int64_t n, int64_t k) {
    return log10Gamma(n+1) - log10Gamma(k+1) - log10Gamma(n-k+1);
}

double log10BinomialProbability(int64_t n, int64_t k, double log10p) {
    if (log10p == 0.0) {
        return -std::numeric_limits<double>::max();
    }
    auto log10MinP = log10(1.0 - pow(10, log10p));
    return log10BinomialCoefficient(n, k) + log10p*static_cast<double>(k) + log10MinP*static_cast<double>(n-k);
}

double log10QualEmpiricalLikelihood(double empiricalQual, double observations, int64_t mismatches) {
    if (observations == 0) {
        return 0.0;
    }
    auto qualToErrorProbLog10 = empiricalQual / -10.0;
    return log10BinomialProbability(observations, mismatches, qualToErrorProbLog10);
}

double qualityToErrorProbability(double phred)  {
    return pow(10, phred/-10);
}

double qualityToProbability(double phred) {
    return 1 - pow(10, phred/-10);
}


double calculateExpectedErrors(int64_t observations, double reportedQuality) {
    return double(observations) * qualityToErrorProbability(reportedQuality);
}


