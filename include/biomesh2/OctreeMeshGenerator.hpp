#pragma once

#include "biomesh2/Octree.hpp"
#include <vector>
#include <array>
#include <unordered_map>

namespace biomesh2 {

/**
 * @brief Hash function for Point3D to enable use in std::unordered_map
 */
struct Point3DHash {
    std::size_t operator()(const Point3D& p) const {
        // Use a simple hash combination
        std::size_t h1 = std::hash<double>{}(p.x);
        std::size_t h2 = std::hash<double>{}(p.y);
        std::size_t h3 = std::hash<double>{}(p.z);
        
        // Combine hashes using a simple method
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

/**
 * @brief Equality function for Point3D to enable use in std::unordered_map
 */
struct Point3DEqual {
    bool operator()(const Point3D& lhs, const Point3D& rhs) const {
        const double epsilon = 1e-12; // Small tolerance for floating point comparison
        return std::abs(lhs.x - rhs.x) < epsilon &&
               std::abs(lhs.y - rhs.y) < epsilon &&
               std::abs(lhs.z - rhs.z) < epsilon;
    }
};

/**
 * @brief Hexahedral mesh data structure
 */
struct HexMesh {
    std::vector<Point3D> nodes;                     // Unique node coordinates
    std::vector<std::array<int, 8>> elements;       // Element connectivity (8 node indices per hex)
    
    /**
     * @brief Get number of nodes in the mesh
     */
    size_t getNodeCount() const { return nodes.size(); }
    
    /**
     * @brief Get number of elements in the mesh
     */
    size_t getElementCount() const { return elements.size(); }
};

/**
 * @brief Octree to hexahedral mesh generator
 */
class OctreeMeshGenerator {
public:
    /**
     * @brief Generate hexahedral mesh from octree leaf nodes
     * @param octree The octree to generate mesh from
     * @return HexMesh containing nodes and element connectivity
     */
    static HexMesh generateHexMesh(const Octree& octree);

private:
    /**
     * @brief Collect all leaf nodes from the octree
     * @param node Current node to examine
     * @param leafNodes Output vector to store leaf nodes
     */
    static void collectLeafNodes(const OctreeNode* node, 
                                std::vector<const OctreeNode*>& leafNodes);
    
    /**
     * @brief Compute the 8 corner nodes for a given octree cell
     * @param node The octree node (leaf cell)
     * @return Array of 8 corner points in standard hexahedral ordering
     */
    static std::array<Point3D, 8> computeCornerNodes(const OctreeNode* node);
    
    /**
     * @brief Assign unique node indices, with deduplication of shared nodes
     * @param cornerNodes Corner nodes from all elements
     * @param nodeMap Map from Point3D to unique index
     * @param uniqueNodes Vector of unique node coordinates
     * @return Vector of element connectivity arrays
     */
    static std::vector<std::array<int, 8>> assignUniqueNodeIndices(
        const std::vector<std::array<Point3D, 8>>& cornerNodes,
        std::unordered_map<Point3D, int, Point3DHash, Point3DEqual>& nodeMap,
        std::vector<Point3D>& uniqueNodes);
};

} // namespace biomesh2