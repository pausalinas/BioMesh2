#include "biomesh2/OctreeMeshGenerator.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace biomesh2;

int main() {
    std::cout << "GiD Mesh Export Demonstration\n";
    std::cout << "=============================\n\n";
    
    try {
        // Example 1: Simple single-element mesh
        std::cout << "1. Creating and exporting single-element mesh...\n";
        Octree simpleOctree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        HexMesh simpleMesh = OctreeMeshGenerator::generateHexMesh(simpleOctree);
        
        std::cout << "   Generated mesh with " << simpleMesh.getNodeCount() 
                  << " nodes and " << simpleMesh.getElementCount() << " elements\n";
        
        OctreeMeshGenerator::exportToGiD(simpleMesh, "simple_mesh.msh");
        std::cout << "   ✓ Exported to simple_mesh.msh\n\n";
        
        // Example 2: Multi-element mesh from subdivided octree
        std::cout << "2. Creating and exporting multi-element mesh...\n";
        Octree complexOctree(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        complexOctree.subdivide(2); // Create 64 leaf cells
        HexMesh complexMesh = OctreeMeshGenerator::generateHexMesh(complexOctree);
        
        std::cout << "   Generated mesh with " << complexMesh.getNodeCount() 
                  << " nodes and " << complexMesh.getElementCount() << " elements\n";
        
        OctreeMeshGenerator::exportToGiD(complexMesh, "complex_mesh.msh");
        std::cout << "   ✓ Exported to complex_mesh.msh\n\n";
        
        // Example 3: Mesh from molecular bounding box (simulated)
        std::cout << "3. Creating mesh for molecular domain...\n";
        // Simulate a molecular bounding box (e.g., from PDB parsing)
        double x_min = -5.0, y_min = -3.0, z_min = -2.0;
        double x_max = 5.0, y_max = 3.0, z_max = 8.0;
        
        Octree molecularOctree(x_min, y_min, z_min, 
                              x_max - x_min, y_max - y_min, z_max - z_min);
        molecularOctree.subdivide(1); // Create 8 leaf cells
        HexMesh molecularMesh = OctreeMeshGenerator::generateHexMesh(molecularOctree);
        
        std::cout << "   Domain: (" << x_min << ", " << y_min << ", " << z_min 
                  << ") to (" << x_max << ", " << y_max << ", " << z_max << ")\n";
        std::cout << "   Generated mesh with " << molecularMesh.getNodeCount() 
                  << " nodes and " << molecularMesh.getElementCount() << " elements\n";
        
        OctreeMeshGenerator::exportToGiD(molecularMesh, "molecular_mesh.msh");
        std::cout << "   ✓ Exported to molecular_mesh.msh\n\n";
        
        // Show sample from the simple mesh file
        std::cout << "4. Sample content from simple_mesh.msh:\n";
        std::ifstream file("simple_mesh.msh");
        if (file.is_open()) {
            std::string line;
            int lineCount = 0;
            while (std::getline(file, line) && lineCount < 15) {
                std::cout << "   " << line << "\n";
                lineCount++;
            }
            if (lineCount >= 15) {
                std::cout << "   ... (truncated)\n";
            }
            file.close();
        }
        
        std::cout << "\nDemo completed successfully!\n";
        std::cout << "Generated mesh files can be imported into GiD or other FEM software.\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}