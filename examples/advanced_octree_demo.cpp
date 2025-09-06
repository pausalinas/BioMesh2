#include "biomesh2/Octree.hpp"
#include "biomesh2/Atom.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>

using namespace biomesh2;

int main() {
    std::cout << "Advanced Octree Features Demonstration" << std::endl;
    std::cout << "======================================" << std::endl;
    
    // Create some sample atoms for testing
    std::vector<std::unique_ptr<Atom>> atoms;
    
    // Add atoms in specific regions to test adaptive refinement
    auto atom1 = std::make_unique<Atom>("C", 0.67, 12.0);
    atom1->setCoordinates(0.25, 0.25, 0.25);
    atoms.push_back(std::move(atom1));
    
    auto atom2 = std::make_unique<Atom>("N", 0.56, 14.0);
    atom2->setCoordinates(0.26, 0.26, 0.26);
    atoms.push_back(std::move(atom2));
    
    auto atom3 = std::make_unique<Atom>("O", 0.48, 16.0);
    atom3->setCoordinates(0.27, 0.27, 0.27);
    atoms.push_back(std::move(atom3));
    
    auto atom4 = std::make_unique<Atom>("H", 0.31, 1.0);
    atom4->setCoordinates(0.75, 0.75, 0.75);
    atoms.push_back(std::move(atom4));
    
    std::cout << "\n1. Creating octree with " << atoms.size() << " atoms..." << std::endl;
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    
    std::cout << "\n2. Demonstrating adaptive subdivision based on atom density..." << std::endl;
    octree.adaptiveSubdivide(atoms, 4, 0.01, RefinementCriterion::ATOM_DENSITY, 2.0, true, true);
    
    auto stats = octree.getTreeStatistics();
    std::cout << "After adaptive subdivision:" << std::endl;
    std::cout << "  - Total nodes: " << static_cast<size_t>(stats["Total Nodes"]) << std::endl;
    std::cout << "  - Leaf nodes: " << static_cast<size_t>(stats["Leaf Nodes"]) << std::endl;
    std::cout << "  - Max depth: " << static_cast<int>(stats["Max Depth"]) << std::endl;
    std::cout << "  - Memory usage: " << static_cast<size_t>(stats["Memory Usage (bytes)"]) << " bytes" << std::endl;
    std::cout << "  - Is 2:1 balanced: " << (stats["Is 2:1 Balanced"] > 0.5 ? "Yes" : "No") << std::endl;
    
    std::cout << "\n3. Showing tree structure with adaptive refinement (depth 0-2):" << std::endl;
    octree.printTree(2, false);
    
    std::cout << "\n4. Demonstrating node lookup by ID..." << std::endl;
    const OctreeNode& root = octree.getRoot();
    size_t rootId = root.nodeId;
    OctreeNode* foundNode = octree.findNodeById(rootId);
    if (foundNode) {
        std::cout << "Found root node by ID " << rootId << ": ";
        std::cout << "Center at (" << foundNode->center.x << ", " 
                  << foundNode->center.y << ", " << foundNode->center.z << ")" << std::endl;
    }
    
    std::cout << "\n5. Testing neighbor finding..." << std::endl;
    if (!root.isLeaf && root.children[0]) {
        auto faceNeighbors = octree.findFaceNeighbors(root.children[0].get());
        std::cout << "Child 0 has neighbors in " << faceNeighbors.size() << " directions" << std::endl;
        
        for (const auto& pair : faceNeighbors) {
            std::cout << "  Direction " << static_cast<int>(pair.first) 
                      << ": " << pair.second.size() << " neighbors" << std::endl;
        }
    }
    
    std::cout << "\n6. Demonstrating validation..." << std::endl;
    auto errors = octree.validateTree();
    if (errors.empty()) {
        std::cout << "Tree validation: PASSED" << std::endl;
    } else {
        std::cout << "Tree validation: FAILED" << std::endl;
        for (const auto& error : errors) {
            std::cout << "  Error: " << error << std::endl;
        }
    }
    
    std::cout << "\n7. Testing FEM-ready leaf collection..." << std::endl;
    auto femLeaves = octree.getFEMReadyLeaves();
    std::cout << "Found " << femLeaves.size() << " FEM-ready leaf nodes" << std::endl;
    
    std::cout << "\n8. Demonstrating boundary conformity..." << std::endl;
    // Show which nodes intersect with molecular boundaries
    size_t boundaryNodes = 0;
    for (auto* leaf : femLeaves) {
        if (leaf->intersectsBoundary) {
            boundaryNodes++;
        }
    }
    std::cout << "Nodes intersecting molecular boundaries: " << boundaryNodes << std::endl;
    
    std::cout << "\n9. Comparing with uniform subdivision..." << std::endl;
    Octree uniformOctree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    uniformOctree.subdivide(3); // Uniform to depth 3
    
    auto uniformStats = uniformOctree.getTreeStatistics();
    std::cout << "Uniform subdivision (depth 3):" << std::endl;
    std::cout << "  - Total nodes: " << static_cast<size_t>(uniformStats["Total Nodes"]) << std::endl;
    std::cout << "  - Leaf nodes: " << static_cast<size_t>(uniformStats["Leaf Nodes"]) << std::endl;
    std::cout << "  - Memory usage: " << static_cast<size_t>(uniformStats["Memory Usage (bytes)"]) << " bytes" << std::endl;
    
    double adaptiveMemory = stats["Memory Usage (bytes)"];
    double uniformMemory = uniformStats["Memory Usage (bytes)"];
    double savings = ((uniformMemory - adaptiveMemory) / uniformMemory) * 100.0;
    
    std::cout << "\nMemory efficiency: ";
    if (savings > 0) {
        std::cout << "Adaptive saves " << std::fixed << std::setprecision(1) 
                  << savings << "% memory compared to uniform" << std::endl;
    } else {
        std::cout << "Adaptive uses " << std::fixed << std::setprecision(1) 
                  << -savings << "% more memory than uniform" << std::endl;
    }
    
    std::cout << "\nDemonstration complete!" << std::endl;
    std::cout << "\nKey features demonstrated:" << std::endl;
    std::cout << "  ✓ Adaptive refinement based on atom density" << std::endl;
    std::cout << "  ✓ 2:1 balance constraint enforcement" << std::endl;
    std::cout << "  ✓ Node ID tracking and lookup" << std::endl;
    std::cout << "  ✓ Neighbor finding algorithms" << std::endl;
    std::cout << "  ✓ Tree structure validation" << std::endl;
    std::cout << "  ✓ FEM-ready leaf identification" << std::endl;
    std::cout << "  ✓ Boundary conformity detection" << std::endl;
    std::cout << "  ✓ Memory usage optimization" << std::endl;
    std::cout << "  ✓ Enhanced statistics and diagnostics" << std::endl;
    
    return 0;
}