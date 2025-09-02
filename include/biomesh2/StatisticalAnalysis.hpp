#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <array>

namespace biomesh2 {

/**
 * @brief Statistical summary for a dataset
 */
struct StatisticalSummary {
    double mean;
    double median;
    double q25;        // First quartile (25th percentile)
    double q75;        // Third quartile (75th percentile)
    double iqr;        // Interquartile range (Q75 - Q25)
    double min;
    double max;
    double stddev;
    size_t count;
    std::vector<double> outliers;  // Values outside 1.5*IQR from quartiles
    
    StatisticalSummary() : mean(0), median(0), q25(0), q75(0), iqr(0), 
                          min(0), max(0), stddev(0), count(0) {}
};

/**
 * @brief Mann-Whitney U test result
 */
struct MannWhitneyResult {
    double uStatistic;     // U statistic
    double zScore;         // Z-score for large sample approximation
    double pValue;         // P-value (two-tailed)
    bool isSignificant;    // Whether difference is significant at α = 0.05
    double effectSize;     // Effect size (r = Z / sqrt(N))
    
    MannWhitneyResult() : uStatistic(0), zScore(0), pValue(1.0), 
                         isSignificant(false), effectSize(0) {}
};

/**
 * @brief Statistical analysis utilities for stress distribution comparison
 */
class StatisticalAnalysis {
public:
    /**
     * @brief Compute comprehensive statistical summary for a dataset
     * @param data Input dataset
     * @return Statistical summary including median, IQR, outliers, etc.
     */
    static StatisticalSummary computeStatisticalSummary(const std::vector<double>& data);
    
    /**
     * @brief Perform Mann-Whitney U test to compare two independent groups
     * @param group1 First group data (e.g., hotspot stress changes)
     * @param group2 Second group data (e.g., random interface stress changes)
     * @return Mann-Whitney U test results including significance
     */
    static MannWhitneyResult mannWhitneyUTest(const std::vector<double>& group1,
                                             const std::vector<double>& group2);
    
    /**
     * @brief Extract stress change values by location type for statistical analysis
     * @param relativeChanges Map of element ID to relative stress change (RJ2v)
     * @param hotspotElements Element IDs corresponding to hotspot locations
     * @param randomElements Element IDs corresponding to random interface locations
     * @param hotspotChanges Output vector for hotspot stress changes
     * @param randomChanges Output vector for random location stress changes
     */
    static void extractStressChangesByLocation(
        const std::unordered_map<size_t, double>& relativeChanges,
        const std::vector<size_t>& hotspotElements,
        const std::vector<size_t>& randomElements,
        std::vector<double>& hotspotChanges,
        std::vector<double>& randomChanges);
    
    /**
     * @brief Generate boxplot data for visualization
     * @param summary Statistical summary
     * @return Boxplot data as [min, Q25, median, Q75, max, outliers...]
     */
    static std::vector<double> generateBoxplotData(const StatisticalSummary& summary);
    
    /**
     * @brief Compute percentile of a dataset
     * @param data Sorted dataset
     * @param percentile Percentile to compute (0-100)
     * @return Value at the specified percentile
     */
    static double computePercentile(const std::vector<double>& data, double percentile);
    
    /**
     * @brief Identify outliers using the 1.5*IQR rule
     * @param data Input dataset
     * @param q25 First quartile
     * @param q75 Third quartile
     * @return Vector of outlier values
     */
    static std::vector<double> identifyOutliers(const std::vector<double>& data,
                                               double q25, double q75);
    
    /**
     * @brief Compare stress distributions between mutation scenarios
     * @param mut1Results Relative stress changes for MUT1 scenario
     * @param mut2Results Relative stress changes for MUT2 scenario
     * @param hotspotElements Element IDs for hotspot locations
     * @param randomElements Element IDs for random interface locations
     * @return Comprehensive comparison results
     */
    struct DistributionComparison {
        StatisticalSummary mut1Hotspots;
        StatisticalSummary mut1Random;
        StatisticalSummary mut2Hotspots;
        StatisticalSummary mut2Random;
        MannWhitneyResult mut1HotspotsVsRandom;
        MannWhitneyResult mut2HotspotsVsRandom;
        MannWhitneyResult hotspotsM1VsM2;
        MannWhitneyResult randomM1VsM2;
    };
    
    static DistributionComparison compareStressDistributions(
        const std::unordered_map<size_t, double>& mut1Results,
        const std::unordered_map<size_t, double>& mut2Results,
        const std::vector<size_t>& hotspotElements,
        const std::vector<size_t>& randomElements);

private:
    /**
     * @brief Compute U statistic for Mann-Whitney test
     * @param ranks1 Ranks for group 1
     * @param n1 Size of group 1
     * @param n2 Size of group 2
     * @return U statistic
     */
    static double computeUStatistic(const std::vector<double>& ranks1, size_t n1, size_t n2);
    
    /**
     * @brief Assign ranks to combined dataset (handling ties)
     * @param combined Combined and sorted dataset
     * @param group1Indices Indices of group 1 elements in combined dataset
     * @param group2Indices Indices of group 2 elements in combined dataset
     * @param ranks1 Output ranks for group 1
     * @param ranks2 Output ranks for group 2
     */
    static void assignRanks(const std::vector<double>& combined,
                           const std::vector<size_t>& group1Indices,
                           const std::vector<size_t>& group2Indices,
                           std::vector<double>& ranks1,
                           std::vector<double>& ranks2);
    
    /**
     * @brief Compute Z-score for large sample Mann-Whitney approximation
     * @param u U statistic
     * @param n1 Size of group 1
     * @param n2 Size of group 2
     * @return Z-score
     */
    static double computeZScore(double u, size_t n1, size_t n2);
    
    /**
     * @brief Convert Z-score to two-tailed p-value using normal approximation
     * @param z Z-score
     * @return Two-tailed p-value
     */
    static double zScoreToTwoTailedP(double z);
    
    /**
     * @brief Standard normal cumulative distribution function approximation
     * @param z Z-score
     * @return Cumulative probability
     */
    static double normalCDF(double z);
};

} // namespace biomesh2