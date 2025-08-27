#include "biomesh2/OctreeMeshGenerator.hpp"
#include <iostream>
#include <cassert>

/**
 * @brief Test function to verify mesh generation correctness
 */
void testMeshGeneration() {
    using namespace biomesh2;
    
    std::cout << "Running mesh generation tests...\n";
    
    // Test 1: Single cell octree
    {
        std::cout << "Test 1: Single cell octree...\n";
        Octree octree(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
        
        assert(mesh.getElementCount() == 1);
        assert(mesh.getNodeCount() == 8);
        
        // Verify corner coordinates
        const auto& element = mesh.elements[0];
        std::vector<Point3D> expectedCorners = {
            Point3D(0.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0), Point3D(2.0, 2.0, 0.0), Point3D(0.0, 2.0, 0.0),
            Point3D(0.0, 0.0, 2.0), Point3D(2.0, 0.0, 2.0), Point3D(2.0, 2.0, 2.0), Point3D(0.0, 2.0, 2.0)
        };
        
        for (size_t i = 0; i < 8; ++i) {
            const Point3D& actual = mesh.nodes[element[i]];
            const Point3D& expected = expectedCorners[i];
            assert(std::abs(actual.x - expected.x) < 1e-10);
            assert(std::abs(actual.y - expected.y) < 1e-10);
            assert(std::abs(actual.z - expected.z) < 1e-10);
        }
        std::cout << "  ✓ Single cell test passed\n";
    }
    
    // Test 2: 2x2x2 subdivision
    {
        std::cout << "Test 2: 2x2x2 subdivision...\n";
        Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        octree.subdivide(1); // Creates 8 leaf cells
        
        HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
        
        assert(mesh.getElementCount() == 8); // 8 leaf cells
        assert(mesh.getNodeCount() == 27);   // 3x3x3 grid of nodes
        
        std::cout << "  ✓ 2x2x2 subdivision test passed\n";
    }
    
    // Test 3: Empty octree (shouldn't crash)
    {
        std::cout << "Test 3: Verify node sharing...\n";
        Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        octree.subdivide(2); // Creates 64 leaf cells
        
        HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
        
        assert(mesh.getElementCount() == 64);  // 4x4x4 leaf cells
        assert(mesh.getNodeCount() == 125);    // 5x5x5 grid of nodes
        
        // Verify that adjacent elements share nodes
        bool foundSharedNode = false;
        for (size_t i = 0; i < mesh.getElementCount() && !foundSharedNode; ++i) {
            for (size_t j = i + 1; j < mesh.getElementCount() && !foundSharedNode; ++j) {
                const auto& elem1 = mesh.elements[i];
                const auto& elem2 = mesh.elements[j];
                
                // Check if elements share any nodes
                for (int node1 : elem1) {
                    for (int node2 : elem2) {
                        if (node1 == node2) {
                            foundSharedNode = true;
                            break;
                        }
                    }
                    if (foundSharedNode) break;
                }
            }
        }
        assert(foundSharedNode); // Adjacent elements should share nodes
        
        std::cout << "  ✓ Node sharing test passed\n";
    }
    
    std::cout << "All tests passed! ✓\n\n";
}

int main() {
    testMeshGeneration();
    return 0;
}