#pragma once

#include "biomesh2/BoundingBox.hpp"
#include <vector>
#include <memory>
#include <functional>
#include <array>

namespace biomesh2 {

/**
 * @brief Octree node representing a subdivision of 3D space
 */
struct OctreeNode {
    Point3D min;         // Minimum corner of the node's bounding box
    Point3D max;         // Maximum corner of the node's bounding box  
    Point3D center;      // Center point of the node
    Point3D halfSize;    // Half-size in each dimension
    int depth;           // Depth level in the tree (root = 0)
    bool isLeaf;         // True if this node has no children
    std::array<std::unique_ptr<OctreeNode>, 8> children; // 8 octant children
    
    /**
     * @brief Constructor for octree node
     * @param minPoint Minimum corner coordinates
     * @param maxPoint Maximum corner coordinates
     * @param nodeDepth Depth level of this node
     */
    OctreeNode(const Point3D& minPoint, const Point3D& maxPoint, int nodeDepth = 0);
    
    /**
     * @brief Check if a point is contained within this node
     * @param point Point to check
     * @return true if point is inside the node's bounds
     */
    bool contains(const Point3D& point) const;
    
    /**
     * @brief Get the volume of this node
     * @return Volume in cubic units
     */
    double getVolume() const;
};

/**
 * @brief Octree data structure for efficient 3D space partitioning
 * 
 * Recursively subdivides 3D space into 8 octants until termination
 * conditions are met (max depth, min size, or custom condition).
 */
class Octree {
public:
    /**
     * @brief Constructor that creates octree from domain bounds
     * @param x_min Minimum X coordinate of domain
     * @param y_min Minimum Y coordinate of domain 
     * @param z_min Minimum Z coordinate of domain
     * @param x_max Maximum X coordinate of domain
     * @param y_max Maximum Y coordinate of domain
     * @param z_max Maximum Z coordinate of domain
     */
    Octree(double x_min, double y_min, double z_min, 
           double x_max, double y_max, double z_max);
    
    /**
     * @brief Constructor that creates octree from existing bounding box
     * @param boundingBox Existing bounding box to use as domain
     */
    explicit Octree(const BoundingBox& boundingBox);
    
    /**
     * @brief Subdivide the octree until termination conditions are met
     * @param maxDepth Maximum depth to subdivide (default: 8)
     * @param minCellSize Minimum cell size before stopping subdivision (default: 0.001)
     * @param occupancyCheck Optional function to check if cell should be subdivided based on content
     */
    void subdivide(int maxDepth = 8, 
                   double minCellSize = 0.001,
                   std::function<bool(const OctreeNode&)> occupancyCheck = nullptr);
    
    /**
     * @brief Get the root node of the octree
     * @return Reference to root node
     */
    const OctreeNode& getRoot() const { return *root_; }
    
    /**
     * @brief Print the tree structure in a readable format
     * @param showDepth Maximum depth to print (default: all)
     */
    void printTree(int showDepth = -1) const;
    
    /**
     * @brief Get total number of nodes in the tree
     * @return Total node count
     */
    size_t getNodeCount() const;
    
    /**
     * @brief Get number of leaf nodes in the tree
     * @return Leaf node count
     */
    size_t getLeafCount() const;
    
    /**
     * @brief Find the leaf node that contains a given point
     * @param point Point to locate
     * @return Pointer to containing leaf node, nullptr if point is outside domain
     */
    const OctreeNode* findLeaf(const Point3D& point) const;

private:
    std::unique_ptr<OctreeNode> root_;
    
    /**
     * @brief Recursive subdivision implementation
     * @param node Node to subdivide
     * @param maxDepth Maximum depth allowed
     * @param minCellSize Minimum cell size allowed
     * @param occupancyCheck Optional occupancy check function
     */
    void subdivideRecursive(OctreeNode* node, 
                           int maxDepth, 
                           double minCellSize,
                           std::function<bool(const OctreeNode&)> occupancyCheck);
    
    /**
     * @brief Print tree structure recursively
     * @param node Current node to print
     * @param prefix String prefix for indentation
     * @param maxDepth Maximum depth to print
     */
    void printTreeRecursive(const OctreeNode* node, const std::string& prefix, int maxDepth) const;
    
    /**
     * @brief Count nodes recursively
     * @param node Current node
     * @return Number of nodes in subtree
     */
    size_t countNodes(const OctreeNode* node) const;
    
    /**
     * @brief Count leaf nodes recursively
     * @param node Current node
     * @return Number of leaf nodes in subtree
     */
    size_t countLeaves(const OctreeNode* node) const;
    
    /**
     * @brief Find leaf containing point recursively
     * @param node Current node to search
     * @param point Point to locate
     * @return Pointer to containing leaf node or nullptr
     */
    const OctreeNode* findLeafRecursive(const OctreeNode* node, const Point3D& point) const;
};

} // namespace biomesh2