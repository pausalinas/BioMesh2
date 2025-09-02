#pragma once

#include "biomesh2/OctreeMeshGenerator.hpp"
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <cmath>

namespace biomesh2 {

/**
 * @brief Von Mises stress tensor representation
 */
struct StressTensor {
    double xx, yy, zz;  // Normal stress components
    double xy, xz, yz;  // Shear stress components
    
    StressTensor() : xx(0), yy(0), zz(0), xy(0), xz(0), yz(0) {}
    StressTensor(double sxx, double syy, double szz, double sxy, double sxz, double syz)
        : xx(sxx), yy(syy), zz(szz), xy(sxy), xz(sxz), yz(syz) {}
};

/**
 * @brief Displacement vector in 3D space
 */
struct Displacement {
    double x, y, z;
    
    Displacement() : x(0), y(0), z(0) {}
    Displacement(double dx, double dy, double dz) : x(dx), y(dy), z(dz) {}
    
    double magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};

/**
 * @brief Finite element analysis results for a single element
 */
struct ElementResult {
    size_t elementId;
    StressTensor stress;
    Displacement displacement;
    double vonMisesStress;
    
    ElementResult() : elementId(0), vonMisesStress(0.0) {}
};

/**
 * @brief Material properties for finite element analysis
 */
struct MaterialProperties {
    double youngModulus;     // Young's modulus (Pa)
    double poissonRatio;     // Poisson's ratio
    double yieldStrength;    // Yield strength (Pa)
    double density;          // Density (kg/m³)
    
    MaterialProperties() 
        : youngModulus(2.1e11), poissonRatio(0.3), yieldStrength(250e6), density(7850.0) {}
    
    MaterialProperties(double E, double nu, double yield, double rho)
        : youngModulus(E), poissonRatio(nu), yieldStrength(yield), density(rho) {}
};

/**
 * @brief Mutation scenario types for comparative analysis
 */
enum class MutationType {
    CONTROL,    // Control/baseline scenario
    MUT1,       // First mutation scenario
    MUT2        // Second mutation scenario
};

/**
 * @brief Analysis location types for statistical comparison
 */
enum class LocationType {
    HOTSPOT,    // Hotspot locations (hs)
    RANDOM      // Randomly selected interface residues (rd)
};

/**
 * @brief Stress analysis results for a complete analysis
 */
struct AnalysisResults {
    std::vector<ElementResult> elementResults;
    std::unordered_map<size_t, double> relativeStressChanges;  // RJ2v values
    MutationType mutationType;
    LocationType locationType;
    MaterialProperties material;
    
    AnalysisResults() : mutationType(MutationType::CONTROL), locationType(LocationType::HOTSPOT) {}
};

/**
 * @brief Finite element stress analysis class
 */
class StressAnalysis {
public:
    /**
     * @brief Constructor with default material properties
     */
    StressAnalysis();
    
    /**
     * @brief Constructor with custom material properties
     * @param material Material properties for analysis
     */
    explicit StressAnalysis(const MaterialProperties& material);
    
    /**
     * @brief Compute von Mises stress from stress tensor
     * @param stress Stress tensor components
     * @return Von Mises equivalent stress
     */
    static double computeVonMisesStress(const StressTensor& stress);
    
    /**
     * @brief Perform stress analysis on hexahedral mesh
     * @param mesh Hexahedral mesh for analysis
     * @param loadConditions Applied loads/boundary conditions
     * @param mutationType Type of mutation scenario
     * @param locationType Type of analysis location
     * @return Complete analysis results
     */
    AnalysisResults performAnalysis(const HexMesh& mesh, 
                                   const std::vector<Point3D>& loadConditions,
                                   MutationType mutationType = MutationType::CONTROL,
                                   LocationType locationType = LocationType::HOTSPOT);
    
    /**
     * @brief Compute relative stress changes between two analyses
     * @param baselineResults Baseline/control analysis results
     * @param mutationResults Mutation scenario analysis results
     * @return Map of element ID to relative change (RJ2v)
     */
    static std::unordered_map<size_t, double> computeRelativeStressChanges(
        const AnalysisResults& baselineResults,
        const AnalysisResults& mutationResults);
    
    /**
     * @brief Assess proximity to structural failure based on yield strength
     * @param results Analysis results to assess
     * @return Vector of safety factors for each element (>1 = safe, <1 = failure risk)
     */
    std::vector<double> assessStructuralFailureRisk(const AnalysisResults& results) const;
    
    /**
     * @brief Set material properties for analysis
     * @param material New material properties
     */
    void setMaterialProperties(const MaterialProperties& material);
    
    /**
     * @brief Get current material properties
     * @return Current material properties
     */
    const MaterialProperties& getMaterialProperties() const;

private:
    MaterialProperties material_;
    
    /**
     * @brief Simplified finite element calculation for demonstration
     * @param mesh Input mesh
     * @param loadConditions Applied loads
     * @return Vector of element results
     */
    std::vector<ElementResult> computeElementStresses(const HexMesh& mesh,
                                                     const std::vector<Point3D>& loadConditions);
    
    /**
     * @brief Compute element stiffness matrix (simplified)
     * @param elementNodes Corner nodes of hexahedral element
     * @return 24x24 element stiffness matrix (simplified representation)
     */
    std::array<std::array<double, 24>, 24> computeElementStiffness(
        const std::array<Point3D, 8>& elementNodes);
};

} // namespace biomesh2