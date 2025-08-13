#include "biomesh2/Octree.hpp"
#include <iostream>
#include <iomanip>

using namespace biomesh2;

int main() {
    std::cout << "Octree Implementation Demonstration" << std::endl;
    std::cout << "===================================" << std::endl;
    
    // Example domain: (0,0,0) to (1,1,1) as specified in requirements
    std::cout << "\n1. Creating octree with domain (0,0,0) to (1,1,1)..." << std::endl;
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    
    // Display root cell information
    const OctreeNode& root = octree.getRoot();
    std::cout << "Root cell created:" << std::endl;
    std::cout << "  - Bounds: (" << std::fixed << std::setprecision(3) 
              << root.min.x << ", " << root.min.y << ", " << root.min.z << ") to ("
              << root.max.x << ", " << root.max.y << ", " << root.max.z << ")" << std::endl;
    std::cout << "  - Center: (" << root.center.x << ", " << root.center.y << ", " << root.center.z << ")" << std::endl;
    std::cout << "  - Half-size: (" << root.halfSize.x << ", " << root.halfSize.y << ", " << root.halfSize.z << ")" << std::endl;
    std::cout << "  - Volume: " << root.getVolume() << " cubic units" << std::endl;
    std::cout << "  - Depth: " << root.depth << std::endl;
    std::cout << "  - Is leaf: " << (root.isLeaf ? "Yes" : "No") << std::endl;
    
    // Subdivide until depth = 3 as specified in requirements
    std::cout << "\n2. Subdividing octree to depth 3..." << std::endl;
    octree.subdivide(3);
    
    std::cout << "Subdivision complete!" << std::endl;
    std::cout << "  - Total nodes: " << octree.getNodeCount() << std::endl;
    std::cout << "  - Leaf nodes: " << octree.getLeafCount() << std::endl;
    
    // Print the tree structure (limited to depth 2 for readability)
    std::cout << "\n3. Tree structure (showing up to depth 2):" << std::endl;
    octree.printTree(2);
    
    // Demonstrate point location functionality
    std::cout << "\n4. Testing point location functionality:" << std::endl;
    std::vector<Point3D> testPoints = {
        Point3D(0.125, 0.125, 0.125),  // Should be in first octant
        Point3D(0.875, 0.875, 0.875),  // Should be in last octant
        Point3D(0.5, 0.5, 0.5),        // Should be in center
        Point3D(-0.1, 0.5, 0.5),       // Outside domain
    };
    
    for (const auto& point : testPoints) {
        std::cout << "Point (" << std::fixed << std::setprecision(3)
                  << point.x << ", " << point.y << ", " << point.z << "): ";
        
        const OctreeNode* leaf = octree.findLeaf(point);
        if (leaf) {
            std::cout << "Found in leaf at depth " << leaf->depth 
                      << " with center (" << leaf->center.x << ", " 
                      << leaf->center.y << ", " << leaf->center.z << ")" << std::endl;
        } else {
            std::cout << "Outside domain" << std::endl;
        }
    }
    
    // Demonstrate termination conditions
    std::cout << "\n5. Demonstrating termination conditions:" << std::endl;
    
    // Example with min cell size
    std::cout << "\n5a. Octree with min cell size = 0.3:" << std::endl;
    Octree octreeMinSize(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    octreeMinSize.subdivide(10, 0.3);  // High max depth, but limited by min size
    std::cout << "  - Total nodes: " << octreeMinSize.getNodeCount() << std::endl;
    std::cout << "  - Leaf nodes: " << octreeMinSize.getLeafCount() << std::endl;
    
    // Example with occupancy check
    std::cout << "\n5b. Octree with occupancy check (only subdivide if contains center point):" << std::endl;
    Octree octreeOccupancy(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    auto occupancyCheck = [](const OctreeNode& node) -> bool {
        return node.contains(Point3D(0.5, 0.5, 0.5));
    };
    octreeOccupancy.subdivide(5, 0.001, occupancyCheck);
    std::cout << "  - Total nodes: " << octreeOccupancy.getNodeCount() << std::endl;
    std::cout << "  - Leaf nodes: " << octreeOccupancy.getLeafCount() << std::endl;
    
    // Example with different domain
    std::cout << "\n6. Octree with different domain (-2,-2,-2) to (2,2,2):" << std::endl;
    Octree octreeLarge(-2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
    octreeLarge.subdivide(2);
    
    const OctreeNode& largeRoot = octreeLarge.getRoot();
    std::cout << "  - Root volume: " << largeRoot.getVolume() << " cubic units" << std::endl;
    std::cout << "  - Total nodes: " << octreeLarge.getNodeCount() << std::endl;
    std::cout << "  - Each depth-2 leaf volume: " << (largeRoot.getVolume() / octreeLarge.getLeafCount()) << " cubic units" << std::endl;
    
    std::cout << "\nDemonstration complete!" << std::endl;
    
    return 0;
}