#include "biomesh2/Octree.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace biomesh2 {

// OctreeNode implementation
OctreeNode::OctreeNode(const Point3D& minPoint, const Point3D& maxPoint, int nodeDepth)
    : min(minPoint), max(maxPoint), depth(nodeDepth), isLeaf(true) {
    
    // Calculate center point
    center = Point3D(
        (min.x + max.x) * 0.5,
        (min.y + max.y) * 0.5,
        (min.z + max.z) * 0.5
    );
    
    // Calculate half-size in each dimension
    halfSize = Point3D(
        (max.x - min.x) * 0.5,
        (max.y - min.y) * 0.5,
        (max.z - min.z) * 0.5
    );
    
    // Initialize all children to nullptr
    for (auto& child : children) {
        child = nullptr;
    }
}

bool OctreeNode::contains(const Point3D& point) const {
    return (point.x >= min.x && point.x <= max.x &&
            point.y >= min.y && point.y <= max.y &&
            point.z >= min.z && point.z <= max.z);
}

double OctreeNode::getVolume() const {
    return (max.x - min.x) * (max.y - min.y) * (max.z - min.z);
}

// Octree implementation
Octree::Octree(double x_min, double y_min, double z_min, 
               double x_max, double y_max, double z_max) {
    Point3D minPoint(x_min, y_min, z_min);
    Point3D maxPoint(x_max, y_max, z_max);
    root_ = std::make_unique<OctreeNode>(minPoint, maxPoint, 0);
}

Octree::Octree(const BoundingBox& boundingBox) {
    Point3D minPoint = boundingBox.getMin();
    Point3D maxPoint = boundingBox.getMax();
    root_ = std::make_unique<OctreeNode>(minPoint, maxPoint, 0);
}

void Octree::subdivide(int maxDepth, double minCellSize, 
                      std::function<bool(const OctreeNode&)> occupancyCheck) {
    if (root_) {
        subdivideRecursive(root_.get(), maxDepth, minCellSize, occupancyCheck);
    }
}

void Octree::subdivideRecursive(OctreeNode* node, int maxDepth, double minCellSize,
                               std::function<bool(const OctreeNode&)> occupancyCheck) {
    if (!node) return;
    
    // Check termination conditions
    if (node->depth >= maxDepth) return;
    
    // Check minimum cell size for potential children
    double minDimension = std::min({
        node->max.x - node->min.x,
        node->max.y - node->min.y,
        node->max.z - node->min.z
    });
    // Don't subdivide if children would be smaller than minCellSize
    if (minDimension * 0.5 < minCellSize) return;
    
    // Check occupancy condition if provided
    if (occupancyCheck && !occupancyCheck(*node)) return;
    
    // Create 8 children (octants)
    double midX = node->center.x;
    double midY = node->center.y;
    double midZ = node->center.z;
    
    // Define the 8 octants in a systematic order
    std::array<std::pair<Point3D, Point3D>, 8> octants = {{
        // Bottom level (z = min to mid)
        {Point3D(node->min.x, node->min.y, node->min.z), Point3D(midX, midY, midZ)},         // 0: ---
        {Point3D(midX, node->min.y, node->min.z), Point3D(node->max.x, midY, midZ)},         // 1: +--
        {Point3D(node->min.x, midY, node->min.z), Point3D(midX, node->max.y, midZ)},         // 2: -+-
        {Point3D(midX, midY, node->min.z), Point3D(node->max.x, node->max.y, midZ)},         // 3: ++-
        
        // Top level (z = mid to max)  
        {Point3D(node->min.x, node->min.y, midZ), Point3D(midX, midY, node->max.z)},         // 4: --+
        {Point3D(midX, node->min.y, midZ), Point3D(node->max.x, midY, node->max.z)},         // 5: +-+
        {Point3D(node->min.x, midY, midZ), Point3D(midX, node->max.y, node->max.z)},         // 6: -++
        {Point3D(midX, midY, midZ), Point3D(node->max.x, node->max.y, node->max.z)}          // 7: +++
    }};
    
    // Create child nodes
    for (size_t i = 0; i < 8; ++i) {
        node->children[i] = std::make_unique<OctreeNode>(
            octants[i].first, octants[i].second, node->depth + 1
        );
    }
    
    // Mark this node as no longer a leaf
    node->isLeaf = false;
    
    // Recursively subdivide children
    for (auto& child : node->children) {
        subdivideRecursive(child.get(), maxDepth, minCellSize, occupancyCheck);
    }
}

void Octree::printTree(int showDepth) const {
    if (!root_) {
        std::cout << "Empty octree" << std::endl;
        return;
    }
    
    std::cout << "Octree Structure:" << std::endl;
    std::cout << "=================" << std::endl;
    printTreeRecursive(root_.get(), "", showDepth);
    
    std::cout << "\nTree Statistics:" << std::endl;
    std::cout << "Total nodes: " << getNodeCount() << std::endl;
    std::cout << "Leaf nodes: " << getLeafCount() << std::endl;
}

void Octree::printTreeRecursive(const OctreeNode* node, const std::string& prefix, int maxDepth) const {
    if (!node || (maxDepth >= 0 && node->depth > maxDepth)) return;
    
    std::cout << prefix;
    std::cout << "├─ Depth " << node->depth;
    std::cout << " | Center: (" << std::fixed << std::setprecision(3) 
              << node->center.x << ", " << node->center.y << ", " << node->center.z << ")";
    std::cout << " | HalfSize: (" << std::fixed << std::setprecision(3)
              << node->halfSize.x << ", " << node->halfSize.y << ", " << node->halfSize.z << ")";
    std::cout << " | Volume: " << std::fixed << std::setprecision(6) << node->getVolume();
    std::cout << " | " << (node->isLeaf ? "LEAF" : "BRANCH");
    std::cout << std::endl;
    
    // Print children if not a leaf
    if (!node->isLeaf) {
        for (size_t i = 0; i < 8; ++i) {
            if (node->children[i]) {
                std::string newPrefix = prefix + "│  ";
                if (i == 7) newPrefix = prefix + "   "; // Last child gets different spacing
                printTreeRecursive(node->children[i].get(), newPrefix, maxDepth);
            }
        }
    }
}

size_t Octree::getNodeCount() const {
    return root_ ? countNodes(root_.get()) : 0;
}

size_t Octree::countNodes(const OctreeNode* node) const {
    if (!node) return 0;
    
    size_t count = 1; // Count this node
    
    // Count all children
    for (const auto& child : node->children) {
        count += countNodes(child.get());
    }
    
    return count;
}

size_t Octree::getLeafCount() const {
    return root_ ? countLeaves(root_.get()) : 0;
}

size_t Octree::countLeaves(const OctreeNode* node) const {
    if (!node) return 0;
    
    if (node->isLeaf) return 1;
    
    size_t count = 0;
    for (const auto& child : node->children) {
        count += countLeaves(child.get());
    }
    
    return count;
}

const OctreeNode* Octree::findLeaf(const Point3D& point) const {
    return root_ ? findLeafRecursive(root_.get(), point) : nullptr;
}

const OctreeNode* Octree::findLeafRecursive(const OctreeNode* node, const Point3D& point) const {
    if (!node || !node->contains(point)) return nullptr;
    
    if (node->isLeaf) return node;
    
    // Search children
    for (const auto& child : node->children) {
        if (auto result = findLeafRecursive(child.get(), point)) {
            return result;
        }
    }
    
    return nullptr;
}

} // namespace biomesh2