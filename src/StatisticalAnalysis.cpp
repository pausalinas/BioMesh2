#include "biomesh2/StatisticalAnalysis.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace biomesh2 {

StatisticalSummary StatisticalAnalysis::computeStatisticalSummary(const std::vector<double>& data) {
    StatisticalSummary summary;
    
    if (data.empty()) {
        return summary;
    }
    
    summary.count = data.size();
    
    // Sort data for percentile calculations
    std::vector<double> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());
    
    summary.min = sortedData.front();
    summary.max = sortedData.back();
    
    // Compute median and quartiles
    summary.median = computePercentile(sortedData, 50.0);
    summary.q25 = computePercentile(sortedData, 25.0);
    summary.q75 = computePercentile(sortedData, 75.0);
    summary.iqr = summary.q75 - summary.q25;
    
    // Compute mean
    summary.mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    
    // Compute standard deviation
    double variance = 0.0;
    for (double value : data) {
        variance += (value - summary.mean) * (value - summary.mean);
    }
    variance /= data.size();
    summary.stddev = std::sqrt(variance);
    
    // Identify outliers using 1.5*IQR rule
    summary.outliers = identifyOutliers(data, summary.q25, summary.q75);
    
    return summary;
}

MannWhitneyResult StatisticalAnalysis::mannWhitneyUTest(const std::vector<double>& group1,
                                                       const std::vector<double>& group2) {
    MannWhitneyResult result;
    
    if (group1.empty() || group2.empty()) {
        return result;
    }
    
    size_t n1 = group1.size();
    size_t n2 = group2.size();
    
    // Combine and sort both groups while tracking group membership
    std::vector<std::pair<double, int>> combined;  // value, group (0 or 1)
    combined.reserve(n1 + n2);
    
    for (double val : group1) {
        combined.emplace_back(val, 0);
    }
    for (double val : group2) {
        combined.emplace_back(val, 1);
    }
    
    std::sort(combined.begin(), combined.end());
    
    // Assign ranks (handling ties with average rank)
    std::vector<double> ranks1, ranks2;
    std::vector<size_t> group1Indices, group2Indices;
    
    for (size_t i = 0; i < combined.size(); ++i) {
        if (combined[i].second == 0) {
            group1Indices.push_back(i);
        } else {
            group2Indices.push_back(i);
        }
    }
    
    // Simplified ranking (for ties, use average rank)
    std::vector<double> allRanks(combined.size());
    for (size_t i = 0; i < combined.size(); ++i) {
        allRanks[i] = i + 1;  // Ranks start from 1
    }
    
    // Handle ties by averaging ranks
    for (size_t i = 0; i < combined.size(); ) {
        size_t j = i;
        while (j < combined.size() && combined[j].first == combined[i].first) {
            ++j;
        }
        
        if (j - i > 1) {  // There are ties
            double avgRank = 0.0;
            for (size_t k = i; k < j; ++k) {
                avgRank += allRanks[k];
            }
            avgRank /= (j - i);
            
            for (size_t k = i; k < j; ++k) {
                allRanks[k] = avgRank;
            }
        }
        i = j;
    }
    
    // Extract ranks for each group
    for (size_t idx : group1Indices) {
        ranks1.push_back(allRanks[idx]);
    }
    for (size_t idx : group2Indices) {
        ranks2.push_back(allRanks[idx]);
    }
    
    // Compute U statistics
    double sumRanks1 = std::accumulate(ranks1.begin(), ranks1.end(), 0.0);
    double u1 = sumRanks1 - (n1 * (n1 + 1)) / 2.0;
    double u2 = n1 * n2 - u1;
    
    result.uStatistic = std::min(u1, u2);
    
    // Compute Z-score for large sample approximation
    result.zScore = computeZScore(result.uStatistic, n1, n2);
    
    // Compute p-value
    result.pValue = zScoreToTwoTailedP(result.zScore);
    result.isSignificant = result.pValue < 0.05;
    
    // Effect size (r = Z / sqrt(N))
    result.effectSize = std::abs(result.zScore) / std::sqrt(n1 + n2);
    
    return result;
}

void StatisticalAnalysis::extractStressChangesByLocation(
    const std::unordered_map<size_t, double>& relativeChanges,
    const std::vector<size_t>& hotspotElements,
    const std::vector<size_t>& randomElements,
    std::vector<double>& hotspotChanges,
    std::vector<double>& randomChanges) {
    
    hotspotChanges.clear();
    randomChanges.clear();
    
    // Extract hotspot stress changes
    for (size_t elementId : hotspotElements) {
        auto it = relativeChanges.find(elementId);
        if (it != relativeChanges.end()) {
            hotspotChanges.push_back(it->second);
        }
    }
    
    // Extract random location stress changes
    for (size_t elementId : randomElements) {
        auto it = relativeChanges.find(elementId);
        if (it != relativeChanges.end()) {
            randomChanges.push_back(it->second);
        }
    }
}

std::vector<double> StatisticalAnalysis::generateBoxplotData(const StatisticalSummary& summary) {
    std::vector<double> boxplotData;
    
    // Standard boxplot data: [min, Q25, median, Q75, max]
    boxplotData.push_back(summary.min);
    boxplotData.push_back(summary.q25);
    boxplotData.push_back(summary.median);
    boxplotData.push_back(summary.q75);
    boxplotData.push_back(summary.max);
    
    // Add outliers
    for (double outlier : summary.outliers) {
        boxplotData.push_back(outlier);
    }
    
    return boxplotData;
}

double StatisticalAnalysis::computePercentile(const std::vector<double>& data, double percentile) {
    if (data.empty()) return 0.0;
    if (data.size() == 1) return data[0];
    
    double index = (percentile / 100.0) * (data.size() - 1);
    size_t lowerIndex = static_cast<size_t>(std::floor(index));
    size_t upperIndex = static_cast<size_t>(std::ceil(index));
    
    if (lowerIndex == upperIndex) {
        return data[lowerIndex];
    }
    
    double fraction = index - lowerIndex;
    return data[lowerIndex] * (1.0 - fraction) + data[upperIndex] * fraction;
}

std::vector<double> StatisticalAnalysis::identifyOutliers(const std::vector<double>& data,
                                                         double q25, double q75) {
    std::vector<double> outliers;
    double iqr = q75 - q25;
    double lowerFence = q25 - 1.5 * iqr;
    double upperFence = q75 + 1.5 * iqr;
    
    for (double value : data) {
        if (value < lowerFence || value > upperFence) {
            outliers.push_back(value);
        }
    }
    
    return outliers;
}

StatisticalAnalysis::DistributionComparison StatisticalAnalysis::compareStressDistributions(
    const std::unordered_map<size_t, double>& mut1Results,
    const std::unordered_map<size_t, double>& mut2Results,
    const std::vector<size_t>& hotspotElements,
    const std::vector<size_t>& randomElements) {
    
    DistributionComparison comparison;
    
    // Extract stress changes by location for MUT1
    std::vector<double> mut1Hotspots, mut1Random;
    extractStressChangesByLocation(mut1Results, hotspotElements, randomElements,
                                 mut1Hotspots, mut1Random);
    
    // Extract stress changes by location for MUT2
    std::vector<double> mut2Hotspots, mut2Random;
    extractStressChangesByLocation(mut2Results, hotspotElements, randomElements,
                                 mut2Hotspots, mut2Random);
    
    // Compute statistical summaries
    comparison.mut1Hotspots = computeStatisticalSummary(mut1Hotspots);
    comparison.mut1Random = computeStatisticalSummary(mut1Random);
    comparison.mut2Hotspots = computeStatisticalSummary(mut2Hotspots);
    comparison.mut2Random = computeStatisticalSummary(mut2Random);
    
    // Perform Mann-Whitney U tests
    comparison.mut1HotspotsVsRandom = mannWhitneyUTest(mut1Hotspots, mut1Random);
    comparison.mut2HotspotsVsRandom = mannWhitneyUTest(mut2Hotspots, mut2Random);
    comparison.hotspotsM1VsM2 = mannWhitneyUTest(mut1Hotspots, mut2Hotspots);
    comparison.randomM1VsM2 = mannWhitneyUTest(mut1Random, mut2Random);
    
    return comparison;
}

double StatisticalAnalysis::computeUStatistic(const std::vector<double>& ranks1, size_t n1, size_t n2) {
    double sumRanks1 = std::accumulate(ranks1.begin(), ranks1.end(), 0.0);
    return sumRanks1 - (n1 * (n1 + 1)) / 2.0;
}

double StatisticalAnalysis::computeZScore(double u, size_t n1, size_t n2) {
    double meanU = (n1 * n2) / 2.0;
    double stdU = std::sqrt((n1 * n2 * (n1 + n2 + 1)) / 12.0);
    
    if (stdU == 0.0) return 0.0;
    
    // Continuity correction for normal approximation
    double correction = (u > meanU) ? -0.5 : 0.5;
    return (u - meanU + correction) / stdU;
}

double StatisticalAnalysis::zScoreToTwoTailedP(double z) {
    return 2.0 * (1.0 - normalCDF(std::abs(z)));
}

double StatisticalAnalysis::normalCDF(double z) {
    // Approximation using the complementary error function
    return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
}

} // namespace biomesh2