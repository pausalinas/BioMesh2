#include "biomesh2/OctreeMeshGenerator.hpp"
#include <cmath>

namespace biomesh2 {

HexMesh OctreeMeshGenerator::generateHexMesh(const Octree& octree) {
    HexMesh mesh;
    
    // Step 1: Collect all leaf nodes
    std::vector<const OctreeNode*> leafNodes;
    collectLeafNodes(&octree.getRoot(), leafNodes);
    
    if (leafNodes.empty()) {
        return mesh; // Return empty mesh if no leaf nodes
    }
    
    // Step 2: Compute corner nodes for all leaf cells
    std::vector<std::array<Point3D, 8>> allCornerNodes;
    allCornerNodes.reserve(leafNodes.size());
    
    for (const OctreeNode* leaf : leafNodes) {
        allCornerNodes.push_back(computeCornerNodes(leaf));
    }
    
    // Step 3: Assign unique node indices with deduplication
    std::unordered_map<Point3D, int, Point3DHash, Point3DEqual> nodeMap;
    mesh.elements = assignUniqueNodeIndices(allCornerNodes, nodeMap, mesh.nodes);
    
    return mesh;
}

void OctreeMeshGenerator::collectLeafNodes(const OctreeNode* node, 
                                          std::vector<const OctreeNode*>& leafNodes) {
    if (!node) return;
    
    if (node->isLeaf) {
        leafNodes.push_back(node);
        return;
    }
    
    // Recursively collect from all children
    for (const auto& child : node->children) {
        if (child) {
            collectLeafNodes(child.get(), leafNodes);
        }
    }
}

std::array<Point3D, 8> OctreeMeshGenerator::computeCornerNodes(const OctreeNode* node) {
    // Standard hexahedral element node ordering:
    // 0: (min.x, min.y, min.z) - bottom-left-back
    // 1: (max.x, min.y, min.z) - bottom-right-back  
    // 2: (max.x, max.y, min.z) - bottom-right-front
    // 3: (min.x, max.y, min.z) - bottom-left-front
    // 4: (min.x, min.y, max.z) - top-left-back
    // 5: (max.x, min.y, max.z) - top-right-back
    // 6: (max.x, max.y, max.z) - top-right-front
    // 7: (min.x, max.y, max.z) - top-left-front
    
    return {{
        Point3D(node->min.x, node->min.y, node->min.z), // 0
        Point3D(node->max.x, node->min.y, node->min.z), // 1
        Point3D(node->max.x, node->max.y, node->min.z), // 2
        Point3D(node->min.x, node->max.y, node->min.z), // 3
        Point3D(node->min.x, node->min.y, node->max.z), // 4
        Point3D(node->max.x, node->min.y, node->max.z), // 5
        Point3D(node->max.x, node->max.y, node->max.z), // 6
        Point3D(node->min.x, node->max.y, node->max.z)  // 7
    }};
}

std::vector<std::array<int, 8>> OctreeMeshGenerator::assignUniqueNodeIndices(
    const std::vector<std::array<Point3D, 8>>& cornerNodes,
    std::unordered_map<Point3D, int, Point3DHash, Point3DEqual>& nodeMap,
    std::vector<Point3D>& uniqueNodes) {
    
    std::vector<std::array<int, 8>> elements;
    elements.reserve(cornerNodes.size());
    
    int nextNodeIndex = 0;
    
    for (const auto& elementCorners : cornerNodes) {
        std::array<int, 8> elementConnectivity;
        
        // For each corner of this element
        for (size_t i = 0; i < 8; ++i) {
            const Point3D& corner = elementCorners[i];
            
            // Check if this node already exists
            auto it = nodeMap.find(corner);
            if (it != nodeMap.end()) {
                // Node already exists, reuse its index
                elementConnectivity[i] = it->second;
            } else {
                // New node, assign it a new index
                elementConnectivity[i] = nextNodeIndex;
                nodeMap[corner] = nextNodeIndex;
                uniqueNodes.push_back(corner);
                ++nextNodeIndex;
            }
        }
        
        elements.push_back(elementConnectivity);
    }
    
    return elements;
}

} // namespace biomesh2