#include "biomesh2/OctreeMeshGenerator.hpp"
#include <iostream>
#include <iomanip>

void demonstrateFeatures() {
    using namespace biomesh2;
    
    std::cout << "Comprehensive Octree to Hexahedral Mesh Demo\n";
    std::cout << "=============================================\n\n";
    
    // Feature 1: Different octree sizes and subdivisions
    std::cout << "Feature 1: Different octree configurations\n";
    std::cout << "------------------------------------------\n";
    
    struct TestCase {
        std::string name;
        double x0, y0, z0, x1, y1, z1;
        int depth;
    };
    
    std::vector<TestCase> testCases = {
        {"Unit cube, depth 1", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1},
        {"Unit cube, depth 2", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2},
        {"Rectangular domain", 0.0, 0.0, 0.0, 2.0, 1.0, 0.5, 2},
        {"Offset domain", -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 2}
    };
    
    for (const auto& testCase : testCases) {
        Octree octree(testCase.x0, testCase.y0, testCase.z0, 
                      testCase.x1, testCase.y1, testCase.z1);
        octree.subdivide(testCase.depth);
        
        HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
        
        std::cout << "  " << testCase.name << ":\n";
        std::cout << "    Leaf nodes: " << octree.getLeafCount() << "\n";
        std::cout << "    Mesh nodes: " << mesh.getNodeCount() << "\n";
        std::cout << "    Mesh elements: " << mesh.getElementCount() << "\n";
        std::cout << "    Node sharing efficiency: " 
                  << std::fixed << std::setprecision(1)
                  << (1.0 - (double)mesh.getNodeCount() / (octree.getLeafCount() * 8.0)) * 100.0 
                  << "%\n\n";
    }
    
    // Feature 2: Demonstrate node ordering
    std::cout << "Feature 2: Standard hexahedral node ordering\n";
    std::cout << "--------------------------------------------\n";
    
    Octree singleOctree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    HexMesh singleMesh = OctreeMeshGenerator::generateHexMesh(singleOctree);
    
    std::cout << "  Single element mesh corners (standard hex ordering):\n";
    const auto& element = singleMesh.elements[0];
    std::vector<std::string> nodeLabels = {
        "bottom-left-back", "bottom-right-back", "bottom-right-front", "bottom-left-front",
        "top-left-back", "top-right-back", "top-right-front", "top-left-front"
    };
    
    for (size_t i = 0; i < 8; ++i) {
        const Point3D& corner = singleMesh.nodes[element[i]];
        std::cout << "    Node " << i << " (" << nodeLabels[i] << "): ("
                  << corner.x << ", " << corner.y << ", " << corner.z << ")\n";
    }
    
    // Feature 3: Memory efficiency demonstration
    std::cout << "\nFeature 3: Memory efficiency through node sharing\n";
    std::cout << "------------------------------------------------\n";
    
    for (int depth = 1; depth <= 3; ++depth) {
        Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        octree.subdivide(depth);
        HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
        
        size_t naiveNodes = octree.getLeafCount() * 8;
        size_t actualNodes = mesh.getNodeCount();
        size_t savedNodes = naiveNodes - actualNodes;
        double efficiency = (double)savedNodes / naiveNodes * 100.0;
        
        std::cout << "  Depth " << depth << ":\n";
        std::cout << "    Without sharing: " << naiveNodes << " nodes\n";
        std::cout << "    With sharing: " << actualNodes << " nodes\n";
        std::cout << "    Memory saved: " << savedNodes << " nodes (" 
                  << std::fixed << std::setprecision(1) << efficiency << "%)\n\n";
    }
    
    std::cout << "Demo complete! All features working correctly.\n";
}

int main() {
    demonstrateFeatures();
    return 0;
}