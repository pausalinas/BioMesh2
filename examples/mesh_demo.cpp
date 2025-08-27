#include "biomesh2/OctreeMeshGenerator.hpp"
#include <iostream>
#include <iomanip>

int main() {
    using namespace biomesh2;
    
    std::cout << "Octree to Hexahedral Mesh Generation Demo\n";
    std::cout << "==========================================\n\n";
    
    // Create a simple octree
    std::cout << "1. Creating octree with domain (0,0,0) to (1,1,1)...\n";
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    
    // Subdivide to create some leaf nodes
    std::cout << "2. Subdividing octree to depth 2...\n";
    octree.subdivide(2); // Create 64 leaf cells
    
    std::cout << "   - Total nodes: " << octree.getNodeCount() << "\n";
    std::cout << "   - Leaf nodes: " << octree.getLeafCount() << "\n\n";
    
    // Generate mesh from octree
    std::cout << "3. Generating hexahedral mesh from octree leaf nodes...\n";
    HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
    
    std::cout << "   Mesh generated successfully!\n";
    std::cout << "   - Number of nodes: " << mesh.getNodeCount() << "\n";
    std::cout << "   - Number of elements: " << mesh.getElementCount() << "\n\n";
    
    // Verify the mesh has the expected number of elements
    size_t expectedElements = octree.getLeafCount();
    if (mesh.getElementCount() == expectedElements) {
        std::cout << "✓ Mesh has correct number of elements (one per leaf node)\n";
    } else {
        std::cout << "✗ Unexpected number of elements! Expected: " << expectedElements 
                  << ", Got: " << mesh.getElementCount() << "\n";
    }
    
    // Show some example nodes and elements
    std::cout << "\n4. Sample mesh data:\n";
    
    // Show first few nodes
    std::cout << "   First 8 nodes:\n";
    for (size_t i = 0; i < std::min(8ul, mesh.getNodeCount()); ++i) {
        const Point3D& node = mesh.nodes[i];
        std::cout << "     Node " << i << ": (" << std::fixed << std::setprecision(3)
                  << node.x << ", " << node.y << ", " << node.z << ")\n";
    }
    
    // Show first few elements
    std::cout << "   First 3 elements:\n";
    for (size_t i = 0; i < std::min(3ul, mesh.getElementCount()); ++i) {
        const auto& element = mesh.elements[i];
        std::cout << "     Element " << i << ": [";
        for (size_t j = 0; j < 8; ++j) {
            std::cout << element[j];
            if (j < 7) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    // Test with a single-cell octree for validation
    std::cout << "\n5. Testing with single-cell octree...\n";
    Octree singleOctree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    // Don't subdivide - should have 1 leaf node
    
    HexMesh singleMesh = OctreeMeshGenerator::generateHexMesh(singleOctree);
    std::cout << "   Single cell mesh:\n";
    std::cout << "   - Nodes: " << singleMesh.getNodeCount() << " (should be 8)\n";
    std::cout << "   - Elements: " << singleMesh.getElementCount() << " (should be 1)\n";
    
    if (singleMesh.getNodeCount() == 8 && singleMesh.getElementCount() == 1) {
        std::cout << "✓ Single cell mesh is correct\n";
        
        // Verify the 8 corners are correct
        std::cout << "   Corner verification:\n";
        const auto& element = singleMesh.elements[0];
        for (size_t i = 0; i < 8; ++i) {
            const Point3D& corner = singleMesh.nodes[element[i]];
            std::cout << "     Corner " << i << ": (" << corner.x << ", " 
                      << corner.y << ", " << corner.z << ")\n";
        }
    } else {
        std::cout << "✗ Single cell mesh has incorrect dimensions\n";
    }
    
    std::cout << "\nDemo complete!\n";
    return 0;
}