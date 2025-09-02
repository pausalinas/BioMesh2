#include "biomesh2/StressAnalysis.hpp"
#include <cmath>
#include <algorithm>
#include <random>

namespace biomesh2 {

StressAnalysis::StressAnalysis() : material_() {
    // Default constructor uses default material properties
}

StressAnalysis::StressAnalysis(const MaterialProperties& material) : material_(material) {
    // Constructor with custom material properties
}

double StressAnalysis::computeVonMisesStress(const StressTensor& stress) {
    // Von Mises stress formula: σvm = √((σxx-σyy)² + (σyy-σzz)² + (σzz-σxx)² + 6(τxy² + τyz² + τzx²))/2
    double term1 = (stress.xx - stress.yy) * (stress.xx - stress.yy);
    double term2 = (stress.yy - stress.zz) * (stress.yy - stress.zz);
    double term3 = (stress.zz - stress.xx) * (stress.zz - stress.xx);
    double term4 = 6.0 * (stress.xy * stress.xy + stress.yz * stress.yz + stress.xz * stress.xz);
    
    return std::sqrt((term1 + term2 + term3 + term4) / 2.0);
}

AnalysisResults StressAnalysis::performAnalysis(const HexMesh& mesh, 
                                               const std::vector<Point3D>& loadConditions,
                                               MutationType mutationType,
                                               LocationType locationType) {
    AnalysisResults results;
    results.mutationType = mutationType;
    results.locationType = locationType;
    results.material = material_;
    
    // Perform finite element analysis
    results.elementResults = computeElementStresses(mesh, loadConditions);
    
    // Compute von Mises stress for each element
    for (auto& elementResult : results.elementResults) {
        elementResult.vonMisesStress = computeVonMisesStress(elementResult.stress);
    }
    
    return results;
}

std::unordered_map<size_t, double> StressAnalysis::computeRelativeStressChanges(
    const AnalysisResults& baselineResults,
    const AnalysisResults& mutationResults) {
    
    std::unordered_map<size_t, double> relativeChanges;
    
    // Create maps for fast lookup
    std::unordered_map<size_t, double> baselineStresses;
    std::unordered_map<size_t, double> mutationStresses;
    
    for (const auto& result : baselineResults.elementResults) {
        baselineStresses[result.elementId] = result.vonMisesStress;
    }
    
    for (const auto& result : mutationResults.elementResults) {
        mutationStresses[result.elementId] = result.vonMisesStress;
    }
    
    // Compute relative changes: RJ2v = (σvm_mutation - σvm_baseline) / σvm_baseline
    for (const auto& baseline : baselineStresses) {
        size_t elementId = baseline.first;
        double baselineStress = baseline.second;
        
        auto mutationIt = mutationStresses.find(elementId);
        if (mutationIt != mutationStresses.end() && baselineStress > 1e-12) {
            double mutationStress = mutationIt->second;
            double relativeChange = (mutationStress - baselineStress) / baselineStress;
            relativeChanges[elementId] = relativeChange;
        }
    }
    
    return relativeChanges;
}

std::vector<double> StressAnalysis::assessStructuralFailureRisk(const AnalysisResults& results) const {
    std::vector<double> safetyFactors;
    safetyFactors.reserve(results.elementResults.size());
    
    for (const auto& elementResult : results.elementResults) {
        // Safety factor = yield strength / von Mises stress
        // Values > 1 indicate safe operation, < 1 indicate failure risk
        double safetyFactor = (elementResult.vonMisesStress > 1e-12) ? 
            material_.yieldStrength / elementResult.vonMisesStress : 
            std::numeric_limits<double>::infinity();
        
        safetyFactors.push_back(safetyFactor);
    }
    
    return safetyFactors;
}

void StressAnalysis::setMaterialProperties(const MaterialProperties& material) {
    material_ = material;
}

const MaterialProperties& StressAnalysis::getMaterialProperties() const {
    return material_;
}

std::vector<ElementResult> StressAnalysis::computeElementStresses(const HexMesh& mesh,
                                                                 const std::vector<Point3D>& loadConditions) {
    std::vector<ElementResult> results;
    results.reserve(mesh.getElementCount());
    
    // This is a simplified finite element analysis for demonstration
    // In a real implementation, this would involve:
    // 1. Assembly of global stiffness matrix
    // 2. Application of boundary conditions
    // 3. Solution of linear system K*u = F
    // 4. Computation of element stresses from displacements
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> stressDist(1e6, 50e6);  // 1-50 MPa range
    std::uniform_real_distribution<> dispDist(-1e-6, 1e-6);  // μm range displacements
    
    for (size_t i = 0; i < mesh.getElementCount(); ++i) {
        ElementResult result;
        result.elementId = i;
        
        // Simplified stress calculation based on element position and applied loads
        const auto& element = mesh.elements[i];
        
        // Compute element centroid
        Point3D centroid(0, 0, 0);
        for (int nodeIdx : element) {
            const Point3D& node = mesh.nodes[nodeIdx];
            centroid.x += node.x;
            centroid.y += node.y;
            centroid.z += node.z;
        }
        centroid.x /= 8.0;
        centroid.y /= 8.0;
        centroid.z /= 8.0;
        
        // Simple stress calculation based on distance to loads and material properties
        double totalLoad = 0.0;
        double avgDistance = 0.0;
        
        for (const auto& load : loadConditions) {
            double dx = centroid.x - load.x;
            double dy = centroid.y - load.y;
            double dz = centroid.z - load.z;
            double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            avgDistance += distance;
            totalLoad += 1.0 / (1.0 + distance);  // Inverse distance weighting
        }
        
        if (!loadConditions.empty()) {
            avgDistance /= loadConditions.size();
            totalLoad /= loadConditions.size();
        }
        
        // Generate realistic stress values based on load distribution
        double baseStress = stressDist(gen) * totalLoad;
        
        // Create stress tensor with some realistic ratios
        result.stress.xx = baseStress;
        result.stress.yy = baseStress * 0.8;
        result.stress.zz = baseStress * 0.6;
        result.stress.xy = baseStress * 0.1;
        result.stress.xz = baseStress * 0.05;
        result.stress.yz = baseStress * 0.05;
        
        // Simplified displacement calculation
        double scaleFactor = baseStress / material_.youngModulus;
        result.displacement.x = dispDist(gen) * scaleFactor;
        result.displacement.y = dispDist(gen) * scaleFactor;
        result.displacement.z = dispDist(gen) * scaleFactor;
        
        results.push_back(result);
    }
    
    return results;
}

std::array<std::array<double, 24>, 24> StressAnalysis::computeElementStiffness(
    const std::array<Point3D, 8>& elementNodes) {
    
    // This is a placeholder for element stiffness matrix computation
    // In a real implementation, this would involve:
    // 1. Numerical integration over element domain
    // 2. Shape function derivatives
    // 3. Jacobian computation and transformation
    // 4. Material property matrix application
    
    std::array<std::array<double, 24>, 24> stiffness = {};
    
    // Initialize with simplified isotropic material stiffness
    double E = material_.youngModulus;
    double nu = material_.poissonRatio;
    double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    
    // This is a highly simplified representation
    // Real implementation would require proper finite element formulation
    for (int i = 0; i < 24; ++i) {
        stiffness[i][i] = factor;  // Diagonal dominance
        for (int j = 0; j < 24; ++j) {
            if (i != j) {
                stiffness[i][j] = factor * nu / (1.0 - nu) * 0.1;  // Off-diagonal coupling
            }
        }
    }
    
    return stiffness;
}

} // namespace biomesh2