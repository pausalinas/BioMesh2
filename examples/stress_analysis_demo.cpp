#include "biomesh2/BioMesh2.hpp"
#include "biomesh2/OctreeMeshGenerator.hpp"
#include "biomesh2/StressAnalysis.hpp"
#include "biomesh2/StatisticalAnalysis.hpp"
#include <iostream>
#include <iomanip>

using namespace biomesh2;

int main() {
    std::cout << "BioMesh2 - Von Mises Stress Analysis Demo" << std::endl;
    std::cout << "=========================================" << std::endl << std::endl;
    
    try {
        // Step 1: Create a molecular structure or use built-in demo data
        std::cout << "Step 1: Setting up molecular structure and mesh generation..." << std::endl;
        
        // Create a simple octree-based mesh for demonstration
        Octree octree(0.0, 0.0, 0.0, 10.0, 10.0, 5.0);  // 10x10x5 Å domain
        octree.subdivide(3);  // Create refined mesh
        
        HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
        std::cout << "   Generated mesh with " << mesh.getNodeCount() << " nodes and " 
                  << mesh.getElementCount() << " elements" << std::endl;
        
        // Step 2: Set up material properties for molecular/protein analysis
        std::cout << "\nStep 2: Configuring material properties..." << std::endl;
        MaterialProperties proteinMaterial(
            5e9,    // Young's modulus: 5 GPa (typical for proteins)
            0.4,    // Poisson's ratio: 0.4 (near incompressible)
            100e6,  // Yield strength: 100 MPa
            1200.0  // Density: 1200 kg/m³ (protein density)
        );
        
        StressAnalysis analyzer(proteinMaterial);
        std::cout << "   Material properties set for protein-like structure" << std::endl;
        
        // Step 3: Define mutation scenarios and loading conditions
        std::cout << "\nStep 3: Setting up mutation scenarios..." << std::endl;
        
        // Define hotspot and random element locations
        std::vector<size_t> hotspotElements = {10, 15, 20, 25, 30};  // Example hotspot locations
        std::vector<size_t> randomElements = {5, 35, 55, 75, 95};    // Random interface residues
        
        // Control scenario - baseline loading
        std::vector<Point3D> controlLoads = {
            Point3D(5.0, 5.0, 5.0)  // Central load
        };
        
        // MUT1 scenario - altered loading pattern
        std::vector<Point3D> mut1Loads = {
            Point3D(3.0, 3.0, 5.0),
            Point3D(7.0, 7.0, 5.0)
        };
        
        // MUT2 scenario - different mutation
        std::vector<Point3D> mut2Loads = {
            Point3D(2.0, 8.0, 5.0),
            Point3D(8.0, 2.0, 5.0)
        };
        
        // Step 4: Perform stress analyses for all scenarios
        std::cout << "\nStep 4: Performing finite element stress analyses..." << std::endl;
        
        auto controlResults = analyzer.performAnalysis(mesh, controlLoads, 
                                                      MutationType::CONTROL, 
                                                      LocationType::HOTSPOT);
        std::cout << "   ✓ Control analysis completed" << std::endl;
        
        auto mut1Results = analyzer.performAnalysis(mesh, mut1Loads, 
                                                   MutationType::MUT1, 
                                                   LocationType::HOTSPOT);
        std::cout << "   ✓ MUT1 analysis completed" << std::endl;
        
        auto mut2Results = analyzer.performAnalysis(mesh, mut2Loads, 
                                                   MutationType::MUT2, 
                                                   LocationType::HOTSPOT);
        std::cout << "   ✓ MUT2 analysis completed" << std::endl;
        
        // Step 5: Compute relative stress changes (RJ2v)
        std::cout << "\nStep 5: Computing relative stress changes..." << std::endl;
        
        auto rj2hs_mut1 = StressAnalysis::computeRelativeStressChanges(controlResults, mut1Results);
        auto rj2rd_mut1 = StressAnalysis::computeRelativeStressChanges(controlResults, mut1Results);
        auto rj2hs_mut2 = StressAnalysis::computeRelativeStressChanges(controlResults, mut2Results);
        auto rj2rd_mut2 = StressAnalysis::computeRelativeStressChanges(controlResults, mut2Results);
        
        std::cout << "   ✓ Relative stress changes computed for both mutations" << std::endl;
        
        // Step 6: Perform statistical analysis using Mann-Whitney U test
        std::cout << "\nStep 6: Statistical analysis of stress distributions..." << std::endl;
        
        auto comparison = StatisticalAnalysis::compareStressDistributions(
            rj2hs_mut1, rj2hs_mut2, hotspotElements, randomElements);
        
        // Step 7: Display results
        std::cout << "\n=== STRESS ANALYSIS RESULTS ===" << std::endl;
        
        // Display von Mises stress statistics
        std::cout << "\nVon Mises Stress Summary:" << std::endl;
        std::cout << "Control scenario - Mean stress: ";
        double controlMeanStress = 0.0;
        for (const auto& result : controlResults.elementResults) {
            controlMeanStress += result.vonMisesStress;
        }
        controlMeanStress /= controlResults.elementResults.size();
        std::cout << std::scientific << std::setprecision(2) << controlMeanStress << " Pa" << std::endl;
        
        // Display relative stress change statistics
        std::cout << "\nMUT1 Hotspots vs Random Interface Residues:" << std::endl;
        std::cout << "   Hotspots - Median: " << std::fixed << std::setprecision(4) 
                  << comparison.mut1Hotspots.median 
                  << ", IQR: " << comparison.mut1Hotspots.iqr << std::endl;
        std::cout << "   Random - Median: " << comparison.mut1Random.median 
                  << ", IQR: " << comparison.mut1Random.iqr << std::endl;
        std::cout << "   Mann-Whitney U test: p = " << std::setprecision(6) 
                  << comparison.mut1HotspotsVsRandom.pValue;
        if (comparison.mut1HotspotsVsRandom.isSignificant) {
            std::cout << " (SIGNIFICANT)";
        } else {
            std::cout << " (not significant)";
        }
        std::cout << std::endl;
        
        std::cout << "\nMUT2 Hotspots vs Random Interface Residues:" << std::endl;
        std::cout << "   Hotspots - Median: " << std::fixed << std::setprecision(4) 
                  << comparison.mut2Hotspots.median 
                  << ", IQR: " << comparison.mut2Hotspots.iqr << std::endl;
        std::cout << "   Random - Median: " << comparison.mut2Random.median 
                  << ", IQR: " << comparison.mut2Random.iqr << std::endl;
        std::cout << "   Mann-Whitney U test: p = " << std::setprecision(6) 
                  << comparison.mut2HotspotsVsRandom.pValue;
        if (comparison.mut2HotspotsVsRandom.isSignificant) {
            std::cout << " (SIGNIFICANT)";
        } else {
            std::cout << " (not significant)";
        }
        std::cout << std::endl;
        
        // Step 8: Structural failure assessment
        std::cout << "\n=== STRUCTURAL FAILURE ASSESSMENT ===" << std::endl;
        
        auto safetyFactors = analyzer.assessStructuralFailureRisk(mut1Results);
        int elementsAtRisk = 0;
        for (double factor : safetyFactors) {
            if (factor < 2.0) {  // Safety factor < 2 indicates potential risk
                elementsAtRisk++;
            }
        }
        
        std::cout << "Elements with safety factor < 2.0: " << elementsAtRisk 
                  << " out of " << safetyFactors.size() << " total elements" << std::endl;
        
        double minSafetyFactor = *std::min_element(safetyFactors.begin(), safetyFactors.end());
        std::cout << "Minimum safety factor: " << std::fixed << std::setprecision(2) 
                  << minSafetyFactor << std::endl;
        
        if (minSafetyFactor < 1.0) {
            std::cout << "⚠️  WARNING: Some elements exceed yield strength!" << std::endl;
        } else {
            std::cout << "✓ All elements are within safe stress limits" << std::endl;
        }
        
        // Step 9: Boxplot data for visualization
        std::cout << "\n=== BOXPLOT DATA FOR VISUALIZATION ===" << std::endl;
        
        auto mut1HotspotBoxplot = StatisticalAnalysis::generateBoxplotData(comparison.mut1Hotspots);
        auto mut1RandomBoxplot = StatisticalAnalysis::generateBoxplotData(comparison.mut1Random);
        
        std::cout << "MUT1 Hotspots boxplot: [min, Q25, median, Q75, max] = [";
        for (size_t i = 0; i < 5 && i < mut1HotspotBoxplot.size(); ++i) {
            std::cout << std::fixed << std::setprecision(4) << mut1HotspotBoxplot[i];
            if (i < 4) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
        
        std::cout << "MUT1 Random boxplot: [min, Q25, median, Q75, max] = [";
        for (size_t i = 0; i < 5 && i < mut1RandomBoxplot.size(); ++i) {
            std::cout << std::fixed << std::setprecision(4) << mut1RandomBoxplot[i];
            if (i < 4) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
        
        std::cout << "\nAnalysis complete! This demonstrates the quantification and" << std::endl;
        std::cout << "visualization of stress distribution across the mesh domain," << std::endl;
        std::cout << "critical for assessing proximity to structural failure." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}