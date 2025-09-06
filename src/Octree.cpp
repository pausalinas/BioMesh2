#include "biomesh2/Octree.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <queue>
#include <set>
#include <cassert>
#include <limits>

namespace biomesh2 {

// Static member initialization
size_t OctreeNode::nextNodeId_ = 0;
std::mutex OctreeNode::nodeIdMutex_;

size_t OctreeNode::generateNodeId() {
    std::lock_guard<std::mutex> lock(nodeIdMutex_);
    return nextNodeId_++;
}

// Enhanced OctreeNode implementation
OctreeNode::OctreeNode(const Point3D& minPoint, const Point3D& maxPoint, 
                       int nodeDepth, OctreeNode* parentNode, int childIdx)
    : min(minPoint), max(maxPoint), depth(nodeDepth), refinementLevel(nodeDepth),
      isLeaf(true), parent(parentNode), childIndex(childIdx),
      refinementCriterion(RefinementCriterion::NONE), refinementValue(0.0),
      needsRefinement(false), needsCoarsening(false), isBeingProcessed(false),
      atomDensity(0.0), intersectsBoundary(false), boundaryDistance(0.0),
      hasViolated2To1Constraint(false) {
    
    // Generate unique node ID
    nodeId = generateNodeId();
    
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
    
    // Initialize neighbor lists
    for (int i = 0; i <= static_cast<int>(NeighborDirection::VERTEX_POS_POS_POS); ++i) {
        neighbors[static_cast<NeighborDirection>(i)] = std::vector<OctreeNode*>();
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

double OctreeNode::getSurfaceArea() const {
    double dx = max.x - min.x;
    double dy = max.y - min.y;
    double dz = max.z - min.z;
    return 2.0 * (dx * dy + dy * dz + dz * dx);
}

bool OctreeNode::intersectsAtom(const Atom& atom) const {
    // Get atom position and radius
    double ax = atom.getX();
    double ay = atom.getY();
    double az = atom.getZ();
    double radius = atom.getAtomicRadius();
    
    // Find closest point on box to atom center
    double closestX = std::max(min.x, std::min(ax, max.x));
    double closestY = std::max(min.y, std::min(ay, max.y));
    double closestZ = std::max(min.z, std::min(az, max.z));
    
    // Calculate distance from atom center to closest point
    double dx = ax - closestX;
    double dy = ay - closestY;
    double dz = az - closestZ;
    double distanceSquared = dx * dx + dy * dy + dz * dz;
    
    // Intersection occurs if distance is less than radius
    return distanceSquared <= (radius * radius);
}

double OctreeNode::calculateAtomDensity(const std::vector<std::unique_ptr<Atom>>& atoms) {
    containedAtomIds.clear();
    size_t atomCount = 0;
    
    for (size_t i = 0; i < atoms.size(); ++i) {
        if (intersectsAtom(*atoms[i])) {
            containedAtomIds.push_back(i);
            atomCount++;
        }
    }
    
    double volume = getVolume();
    atomDensity = (volume > 0.0) ? static_cast<double>(atomCount) / volume : 0.0;
    return atomDensity;
}

void OctreeNode::updateContainedAtoms(const std::vector<std::unique_ptr<Atom>>& atoms) {
    calculateAtomDensity(atoms);
    
    // Update boundary intersection status
    intersectsBoundary = false;
    boundaryDistance = std::numeric_limits<double>::max();
    
    for (const auto& atom : atoms) {
        double distance = std::sqrt(
            std::pow(atom->getX() - center.x, 2) +
            std::pow(atom->getY() - center.y, 2) +
            std::pow(atom->getZ() - center.z, 2)
        ) - atom->getAtomicRadius();
        
        if (distance < getMinEdgeLength() * 0.5) {
            intersectsBoundary = true;
        }
        
        boundaryDistance = std::min(boundaryDistance, distance);
    }
}

bool OctreeNode::shouldRefine(RefinementCriterion criterion, double threshold,
                             const std::vector<std::unique_ptr<Atom>>& atoms) const {
    switch (criterion) {
        case RefinementCriterion::ATOM_DENSITY:
            return atomDensity > threshold;
            
        case RefinementCriterion::BOUNDARY_PROXIMITY:
            return intersectsBoundary || (boundaryDistance < threshold);
            
        case RefinementCriterion::GEOMETRIC_FEATURE: {
            // Refine if node has significant geometric complexity
            double aspectRatio = getMaxEdgeLength() / getMinEdgeLength();
            return aspectRatio > threshold;
        }
        
        case RefinementCriterion::USER_DEFINED:
            return refinementValue > threshold;
            
        case RefinementCriterion::ERROR_INDICATOR:
            return refinementValue > threshold;
            
        default:
            return false;
    }
}

void OctreeNode::addNeighbor(NeighborDirection direction, OctreeNode* neighbor) {
    if (neighbor) {
        neighbors[direction].push_back(neighbor);
    }
}

const std::vector<OctreeNode*>& OctreeNode::getNeighbors(NeighborDirection direction) const {
    auto it = neighbors.find(direction);
    static const std::vector<OctreeNode*> empty;
    return (it != neighbors.end()) ? it->second : empty;
}

bool OctreeNode::violates2To1Constraint() const {
    // Check face neighbors for level difference > 1
    for (int i = 0; i <= 5; ++i) { // Face neighbors only
        NeighborDirection dir = static_cast<NeighborDirection>(i);
        const auto& neighs = getNeighbors(dir);
        
        for (const auto* neighbor : neighs) {
            if (neighbor && std::abs(neighbor->refinementLevel - refinementLevel) > 1) {
                return true;
            }
        }
    }
    return false;
}

double OctreeNode::getMinEdgeLength() const {
    return std::min({max.x - min.x, max.y - min.y, max.z - min.z});
}

double OctreeNode::getMaxEdgeLength() const {
    return std::max({max.x - min.x, max.y - min.y, max.z - min.z});
}

void OctreeNode::setProcessing(bool processing) {
    isBeingProcessed.store(processing);
}

bool OctreeNode::getProcessing() const {
    return isBeingProcessed.load();
}

// Enhanced Octree implementation
Octree::Octree(double x_min, double y_min, double z_min, 
               double x_max, double y_max, double z_max)
    : enableThreadSafety_(true), enableMemoryOptimization_(true), 
      defaultRefinementThreshold_(1.0) {
    Point3D minPoint(x_min, y_min, z_min);
    Point3D maxPoint(x_max, y_max, z_max);
    root_ = std::make_unique<OctreeNode>(minPoint, maxPoint, 0);
    registerNode(root_.get());
}

Octree::Octree(const BoundingBox& boundingBox)
    : enableThreadSafety_(true), enableMemoryOptimization_(true),
      defaultRefinementThreshold_(1.0) {
    Point3D minPoint = boundingBox.getMin();
    Point3D maxPoint = boundingBox.getMax();
    root_ = std::make_unique<OctreeNode>(minPoint, maxPoint, 0);
    registerNode(root_.get());
}

void Octree::subdivide(int maxDepth, double minCellSize, 
                      std::function<bool(const OctreeNode&)> occupancyCheck) {
    if (root_) {
        subdivideRecursive(root_.get(), maxDepth, minCellSize, occupancyCheck);
        updateNeighborConnectivity();
    }
}

void Octree::adaptiveSubdivide(const std::vector<std::unique_ptr<Atom>>& atoms,
                              int maxDepth, double minCellSize,
                              RefinementCriterion criterion, double threshold,
                              bool enforce2To1, bool enableBoundaryRefinement) {
    if (!root_) return;
    
    // Update atom information for all nodes
    updateNodeAtomInfo(root_.get(), atoms);
    
    // Perform adaptive subdivision
    adaptiveSubdivideRecursive(root_.get(), atoms, maxDepth, minCellSize, criterion, threshold);
    
    // Update neighbor connectivity
    updateNeighborConnectivity();
    
    // Enforce 2:1 balance if requested
    if (enforce2To1) {
        enforce2To1Balance(atoms);
    }
    
    // Boundary refinement if requested
    if (enableBoundaryRefinement) {
        performBoundaryRefinement(atoms, maxDepth, minCellSize);
    }
}

void Octree::enforce2To1Balance(const std::vector<std::unique_ptr<Atom>>& atoms) {
    bool constraintViolated = true;
    int iterations = 0;
    const int maxIterations = 10; // Prevent infinite loops
    
    while (constraintViolated && iterations < maxIterations) {
        constraintViolated = false;
        iterations++;
        
        // Find all nodes that violate 2:1 constraint
        std::vector<OctreeNode*> violatingNodes;
        collectViolatingNodes(root_.get(), violatingNodes);
        
        // Refine violating nodes
        for (auto* node : violatingNodes) {
            if (constrain2To1ForNode(node, atoms)) {
                constraintViolated = true;
            }
        }
        
        // Update neighbor connectivity after refinements
        if (constraintViolated) {
            updateNeighborConnectivity();
        }
    }
}

std::unordered_map<NeighborDirection, std::vector<OctreeNode*>> 
Octree::findFaceNeighbors(OctreeNode* node) const {
    std::unordered_map<NeighborDirection, std::vector<OctreeNode*>> faceNeighbors;
    
    if (!node) return faceNeighbors;
    
    // Implementation of neighbor finding algorithm
    // This is a simplified version - full implementation would be more complex
    double nodeSize = node->getMinEdgeLength();
    
    // Check all 6 face directions
    std::array<std::pair<Point3D, NeighborDirection>, 6> faceDirections = {{
        {Point3D(-nodeSize, 0, 0), NeighborDirection::FACE_X_NEG},
        {Point3D(nodeSize, 0, 0), NeighborDirection::FACE_X_POS},
        {Point3D(0, -nodeSize, 0), NeighborDirection::FACE_Y_NEG},
        {Point3D(0, nodeSize, 0), NeighborDirection::FACE_Y_POS},
        {Point3D(0, 0, -nodeSize), NeighborDirection::FACE_Z_NEG},
        {Point3D(0, 0, nodeSize), NeighborDirection::FACE_Z_POS}
    }};
    
    for (const auto& dirPair : faceDirections) {
        Point3D searchPoint = Point3D(
            node->center.x + dirPair.first.x,
            node->center.y + dirPair.first.y,
            node->center.z + dirPair.first.z
        );
        
        const OctreeNode* neighbor = findLeaf(searchPoint);
        if (neighbor && neighbor != node) {
            faceNeighbors[dirPair.second].push_back(const_cast<OctreeNode*>(neighbor));
        }
    }
    
    return faceNeighbors;
}

bool Octree::refineNode(OctreeNode* node, 
                       const std::vector<std::unique_ptr<Atom>>& atoms,
                       RefinementCriterion criterion, double threshold) {
    if (!node || !node->isLeaf) return false;
    
    // Check if node should be refined
    if (!node->shouldRefine(criterion, threshold, atoms)) {
        return false;
    }
    
    // Create children
    createAllChildren(node);
    
    // Update atom information for children
    for (auto& child : node->children) {
        if (child) {
            child->updateContainedAtoms(atoms);
        }
    }
    
    // Update neighbor connectivity
    updateNeighborsAfterSubdivision(node);
    
    return true;
}

void Octree::updateNeighborConnectivity() {
    if (!root_) return;
    
    // Clear existing neighbor information
    clearNeighborConnectivity(root_.get());
    
    // Rebuild neighbor connectivity
    buildNeighborConnectivity(root_.get());
}

std::vector<OctreeNode*> Octree::getNodesAtLevel(int level) const {
    std::vector<OctreeNode*> nodes;
    if (root_) {
        collectNodesAtLevel(root_.get(), level, nodes);
    }
    return nodes;
}

std::vector<OctreeNode*> Octree::getFEMReadyLeaves() const {
    std::vector<OctreeNode*> femLeaves;
    if (root_) {
        collectFEMReadyLeaves(root_.get(), femLeaves);
    }
    return femLeaves;
}

std::vector<std::string> Octree::validateTree() const {
    std::vector<std::string> errors;
    
    if (!root_) {
        errors.push_back("Tree has no root node");
        return errors;
    }
    
    // Validate tree structure
    validateTreeStructure(root_.get(), errors);
    
    // Validate 2:1 constraint
    if (!is2To1Balanced()) {
        errors.push_back("Tree violates 2:1 balance constraint");
    }
    
    // Validate neighbor connectivity
    validateNeighborConnectivity(root_.get(), errors);
    
    return errors;
}

OctreeNode* Octree::findNodeById(size_t nodeId) const {
    auto it = nodeRegistry_.find(nodeId);
    return (it != nodeRegistry_.end()) ? it->second : nullptr;
}

bool Octree::is2To1Balanced() const {
    if (!root_) return true;
    return checkBalanceRecursive(root_.get());
}

size_t Octree::getMemoryUsage() const {
    size_t usage = sizeof(Octree);
    if (root_) {
        usage += calculateMemoryUsage(root_.get());
    }
    usage += nodeRegistry_.size() * (sizeof(size_t) + sizeof(OctreeNode*));
    return usage;
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
    
    // Create all 8 children
    createAllChildren(node);
    
    // Recursively subdivide children
    for (auto& child : node->children) {
        if (child) {
            subdivideRecursive(child.get(), maxDepth, minCellSize, occupancyCheck);
        }
    }
}

void Octree::adaptiveSubdivideRecursive(OctreeNode* node,
                                       const std::vector<std::unique_ptr<Atom>>& atoms,
                                       int maxDepth, double minCellSize,
                                       RefinementCriterion criterion, double threshold) {
    if (!node) return;
    
    // Check termination conditions
    if (node->depth >= maxDepth) return;
    
    // Check minimum cell size
    if (node->getMinEdgeLength() * 0.5 < minCellSize) return;
    
    // Update node information
    node->updateContainedAtoms(atoms);
    
    // Check if node should be refined
    if (!node->shouldRefine(criterion, threshold, atoms)) return;
    
    // Create children
    createAllChildren(node);
    
    // Recursively subdivide children
    for (auto& child : node->children) {
        if (child) {
            adaptiveSubdivideRecursive(child.get(), atoms, maxDepth, minCellSize, criterion, threshold);
        }
    }
}

void Octree::createAllChildren(OctreeNode* parent) {
    if (!parent || !parent->isLeaf) return;
    
    // Create 8 children (octants)
    double midX = parent->center.x;
    double midY = parent->center.y;
    double midZ = parent->center.z;
    
    // Define the 8 octants in a systematic order
    std::array<std::pair<Point3D, Point3D>, 8> octants = {{
        // Bottom level (z = min to mid)
        {Point3D(parent->min.x, parent->min.y, parent->min.z), Point3D(midX, midY, midZ)},         // 0: ---
        {Point3D(midX, parent->min.y, parent->min.z), Point3D(parent->max.x, midY, midZ)},         // 1: +--
        {Point3D(parent->min.x, midY, parent->min.z), Point3D(midX, parent->max.y, midZ)},         // 2: -+-
        {Point3D(midX, midY, parent->min.z), Point3D(parent->max.x, parent->max.y, midZ)},         // 3: ++-
        
        // Top level (z = mid to max)  
        {Point3D(parent->min.x, parent->min.y, midZ), Point3D(midX, midY, parent->max.z)},         // 4: --+
        {Point3D(midX, parent->min.y, midZ), Point3D(parent->max.x, midY, parent->max.z)},         // 5: +-+
        {Point3D(parent->min.x, midY, midZ), Point3D(midX, parent->max.y, parent->max.z)},         // 6: -++
        {Point3D(midX, midY, midZ), Point3D(parent->max.x, parent->max.y, parent->max.z)}          // 7: +++
    }};
    
    // Create child nodes
    for (size_t i = 0; i < 8; ++i) {
        parent->children[i] = std::make_unique<OctreeNode>(
            octants[i].first, octants[i].second, parent->depth + 1, parent, static_cast<int>(i)
        );
        registerNode(parent->children[i].get());
    }
    
    // Mark this node as no longer a leaf
    parent->isLeaf = false;
}

// Helper methods for the new functionality
void Octree::updateNodeAtomInfo(OctreeNode* node, const std::vector<std::unique_ptr<Atom>>& atoms) {
    if (!node) return;
    
    node->updateContainedAtoms(atoms);
    
    for (auto& child : node->children) {
        if (child) {
            updateNodeAtomInfo(child.get(), atoms);
        }
    }
}

void Octree::performBoundaryRefinement(const std::vector<std::unique_ptr<Atom>>& atoms, 
                                      int maxDepth, double minCellSize) {
    // Collect all leaf nodes that intersect boundaries
    std::vector<OctreeNode*> boundaryNodes;
    collectBoundaryNodes(root_.get(), boundaryNodes);
    
    for (auto* node : boundaryNodes) {
        if (node->depth < maxDepth && node->getMinEdgeLength() * 0.5 >= minCellSize) {
            refineNode(node, atoms, RefinementCriterion::BOUNDARY_PROXIMITY, 0.0);
        }
    }
}

void Octree::collectViolatingNodes(OctreeNode* node, std::vector<OctreeNode*>& violatingNodes) {
    if (!node) return;
    
    if (node->violates2To1Constraint()) {
        violatingNodes.push_back(node);
    }
    
    for (auto& child : node->children) {
        if (child) {
            collectViolatingNodes(child.get(), violatingNodes);
        }
    }
}

void Octree::collectBoundaryNodes(OctreeNode* node, std::vector<OctreeNode*>& boundaryNodes) {
    if (!node) return;
    
    if (node->isLeaf && node->intersectsBoundary) {
        boundaryNodes.push_back(node);
    }
    
    for (auto& child : node->children) {
        if (child) {
            collectBoundaryNodes(child.get(), boundaryNodes);
        }
    }
}

bool Octree::constrain2To1ForNode(OctreeNode* node, 
                                 const std::vector<std::unique_ptr<Atom>>& atoms) {
    if (!node || !node->violates2To1Constraint()) return false;
    
    // Find the neighbor with highest refinement level
    int maxNeighborLevel = node->refinementLevel;
    for (int i = 0; i <= 5; ++i) { // Face neighbors only
        NeighborDirection dir = static_cast<NeighborDirection>(i);
        const auto& neighs = node->getNeighbors(dir);
        
        for (const auto* neighbor : neighs) {
            if (neighbor) {
                maxNeighborLevel = std::max(maxNeighborLevel, neighbor->refinementLevel);
            }
        }
    }
    
    // Refine this node if necessary
    if (maxNeighborLevel - node->refinementLevel > 1) {
        return refineNode(node, atoms, RefinementCriterion::USER_DEFINED, 0.0);
    }
    
    return false;
}

void Octree::clearNeighborConnectivity(OctreeNode* node) {
    if (!node) return;
    
    for (auto& pair : node->neighbors) {
        pair.second.clear();
    }
    
    for (auto& child : node->children) {
        if (child) {
            clearNeighborConnectivity(child.get());
        }
    }
}

void Octree::buildNeighborConnectivity(OctreeNode* node) {
    if (!node) return;
    
    // This is a simplified neighbor finding - a full implementation would be more sophisticated
    if (node->isLeaf) {
        auto faceNeighbors = findFaceNeighbors(node);
        for (const auto& pair : faceNeighbors) {
            for (auto* neighbor : pair.second) {
                node->addNeighbor(pair.first, neighbor);
            }
        }
    }
    
    for (auto& child : node->children) {
        if (child) {
            buildNeighborConnectivity(child.get());
        }
    }
}

void Octree::collectNodesAtLevel(OctreeNode* node, int level, std::vector<OctreeNode*>& nodes) const {
    if (!node) return;
    
    if (node->refinementLevel == level) {
        nodes.push_back(node);
    }
    
    for (auto& child : node->children) {
        if (child) {
            collectNodesAtLevel(child.get(), level, nodes);
        }
    }
}

void Octree::collectFEMReadyLeaves(OctreeNode* node, std::vector<OctreeNode*>& femLeaves) const {
    if (!node) return;
    
    if (node->isLeaf && !node->violates2To1Constraint()) {
        femLeaves.push_back(node);
    }
    
    for (auto& child : node->children) {
        if (child) {
            collectFEMReadyLeaves(child.get(), femLeaves);
        }
    }
}

void Octree::validateTreeStructure(OctreeNode* node, std::vector<std::string>& errors) const {
    if (!node) return;
    
    // Check parent-child consistency
    for (size_t i = 0; i < node->children.size(); ++i) {
        if (node->children[i]) {
            if (node->children[i]->parent != node) {
                errors.push_back("Parent-child relationship inconsistency detected");
            }
            if (node->children[i]->childIndex != static_cast<int>(i)) {
                errors.push_back("Child index inconsistency detected");
            }
        }
    }
    
    // Recursively validate children
    for (auto& child : node->children) {
        if (child) {
            validateTreeStructure(child.get(), errors);
        }
    }
}

void Octree::validateNeighborConnectivity(OctreeNode* node, std::vector<std::string>& errors) const {
    if (!node) return;
    
    // Check neighbor reciprocity (simplified check)
    for (const auto& pair : node->neighbors) {
        for (const auto* neighbor : pair.second) {
            if (neighbor) {
                // In a full implementation, would check if this node is in neighbor's neighbor list
            }
        }
    }
    
    for (auto& child : node->children) {
        if (child) {
            validateNeighborConnectivity(child.get(), errors);
        }
    }
}

bool Octree::checkBalanceRecursive(OctreeNode* node) const {
    if (!node) return true;
    
    if (node->violates2To1Constraint()) {
        return false;
    }
    
    for (auto& child : node->children) {
        if (child && !checkBalanceRecursive(child.get())) {
            return false;
        }
    }
    
    return true;
}

size_t Octree::calculateMemoryUsage(OctreeNode* node) const {
    if (!node) return 0;
    
    size_t usage = sizeof(OctreeNode);
    usage += node->containedAtomIds.size() * sizeof(size_t);
    
    for (auto& child : node->children) {
        if (child) {
            usage += calculateMemoryUsage(child.get());
        }
    }
    
    return usage;
}

void Octree::registerNode(OctreeNode* node) {
    if (node) {
        nodeRegistry_[node->nodeId] = node;
    }
}

void Octree::unregisterNode(OctreeNode* node) {
    if (node) {
        nodeRegistry_.erase(node->nodeId);
    }
}

void Octree::printTree(int showDepth, bool showNeighbors) const {
    if (!root_) {
        std::cout << "Empty octree" << std::endl;
        return;
    }
    
    std::cout << "Advanced Octree Structure:" << std::endl;
    std::cout << "=========================" << std::endl;
    printTreeRecursive(root_.get(), "", showDepth, showNeighbors);
    
    auto stats = getTreeStatistics();
    std::cout << "\nEnhanced Tree Statistics:" << std::endl;
    for (const auto& pair : stats) {
        std::cout << pair.first << ": " << pair.second << std::endl;
    }
}

void Octree::printTreeRecursive(const OctreeNode* node, const std::string& prefix, 
                               int maxDepth, bool showNeighbors) const {
    if (!node || (maxDepth >= 0 && node->depth > maxDepth)) return;
    
    std::cout << prefix;
    std::cout << "├─ Depth " << node->depth << " (Ref " << node->refinementLevel << ")";
    std::cout << " | ID: " << node->nodeId;
    std::cout << " | Center: (" << std::fixed << std::setprecision(3) 
              << node->center.x << ", " << node->center.y << ", " << node->center.z << ")";
    std::cout << " | HalfSize: (" << std::fixed << std::setprecision(3)
              << node->halfSize.x << ", " << node->halfSize.y << ", " << node->halfSize.z << ")";
    std::cout << " | Volume: " << std::fixed << std::setprecision(6) << node->getVolume();
    std::cout << " | Atoms: " << node->containedAtomIds.size();
    std::cout << " | Density: " << std::fixed << std::setprecision(3) << node->atomDensity;
    std::cout << " | " << (node->isLeaf ? "LEAF" : "BRANCH");
    
    if (node->intersectsBoundary) {
        std::cout << " | BOUNDARY";
    }
    
    if (node->violates2To1Constraint()) {
        std::cout << " | VIOLATES_2:1";
    }
    
    std::cout << std::endl;
    
    // Show neighbor information if requested
    if (showNeighbors && node->isLeaf) {
        for (int i = 0; i <= 5; ++i) { // Face neighbors only for brevity
            NeighborDirection dir = static_cast<NeighborDirection>(i);
            const auto& neighs = node->getNeighbors(dir);
            if (!neighs.empty()) {
                std::cout << prefix << "│    Neighbors " << i << ": " << neighs.size() << std::endl;
            }
        }
    }
    
    // Print children if not a leaf
    if (!node->isLeaf) {
        for (size_t i = 0; i < 8; ++i) {
            if (node->children[i]) {
                std::string newPrefix = prefix + "│  ";
                if (i == 7) newPrefix = prefix + "   "; // Last child gets different spacing
                printTreeRecursive(node->children[i].get(), newPrefix, maxDepth, showNeighbors);
            }
        }
    }
}

std::vector<size_t> Octree::getNodeCountByDepth() const {
    std::vector<size_t> counts;
    if (root_) {
        countNodesByDepth(root_.get(), counts);
    }
    return counts;
}

std::unordered_map<std::string, double> Octree::getTreeStatistics() const {
    std::unordered_map<std::string, double> stats;
    
    stats["Total Nodes"] = static_cast<double>(getNodeCount());
    stats["Leaf Nodes"] = static_cast<double>(getLeafCount());
    stats["Memory Usage (bytes)"] = static_cast<double>(getMemoryUsage());
    stats["Is 2:1 Balanced"] = is2To1Balanced() ? 1.0 : 0.0;
    
    if (root_) {
        stats["Max Depth"] = static_cast<double>(getMaxDepth(root_.get()));
        stats["Min Leaf Size"] = getMinLeafSize(root_.get());
        stats["Max Leaf Size"] = getMaxLeafSize(root_.get());
        stats["Average Atom Density"] = getAverageAtomDensity(root_.get());
    }
    
    return stats;
}

// Helper methods for statistics
int Octree::getMaxDepth(OctreeNode* node) const {
    if (!node) return -1;
    
    int maxDepth = node->depth;
    for (auto& child : node->children) {
        if (child) {
            maxDepth = std::max(maxDepth, getMaxDepth(child.get()));
        }
    }
    return maxDepth;
}

double Octree::getMinLeafSize(OctreeNode* node) const {
    if (!node) return std::numeric_limits<double>::max();
    
    double minSize = std::numeric_limits<double>::max();
    
    if (node->isLeaf) {
        minSize = node->getMinEdgeLength();
    }
    
    for (auto& child : node->children) {
        if (child) {
            minSize = std::min(minSize, getMinLeafSize(child.get()));
        }
    }
    
    return minSize;
}

double Octree::getMaxLeafSize(OctreeNode* node) const {
    if (!node) return 0.0;
    
    double maxSize = 0.0;
    
    if (node->isLeaf) {
        maxSize = node->getMaxEdgeLength();
    }
    
    for (auto& child : node->children) {
        if (child) {
            maxSize = std::max(maxSize, getMaxLeafSize(child.get()));
        }
    }
    
    return maxSize;
}

double Octree::getAverageAtomDensity(OctreeNode* node) const {
    if (!node) return 0.0;
    
    double totalDensity = 0.0;
    size_t leafCount = 0;
    
    calculateAverageAtomDensityRecursive(node, totalDensity, leafCount);
    
    return leafCount > 0 ? totalDensity / static_cast<double>(leafCount) : 0.0;
}

void Octree::calculateAverageAtomDensityRecursive(OctreeNode* node, double& totalDensity, size_t& leafCount) const {
    if (!node) return;
    
    if (node->isLeaf) {
        totalDensity += node->atomDensity;
        leafCount++;
    }
    
    for (auto& child : node->children) {
        if (child) {
            calculateAverageAtomDensityRecursive(child.get(), totalDensity, leafCount);
        }
    }
}

void Octree::countNodesByDepth(const OctreeNode* node, std::vector<size_t>& counts) const {
    if (!node) return;
    
    // Ensure vector is large enough
    if (static_cast<size_t>(node->depth) >= counts.size()) {
        counts.resize(node->depth + 1, 0);
    }
    
    counts[node->depth]++;
    
    for (const auto& child : node->children) {
        if (child) {
            countNodesByDepth(child.get(), counts);
        }
    }
}

void Octree::updateNeighborsAfterSubdivision(OctreeNode* parent) {
    // This is a simplified implementation
    // Full implementation would be more sophisticated
    if (!parent) return;
    
    // Clear parent's neighbor information since it's no longer a leaf
    for (auto& pair : parent->neighbors) {
        pair.second.clear();
    }
    
    // Update connectivity will be called globally after subdivision
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