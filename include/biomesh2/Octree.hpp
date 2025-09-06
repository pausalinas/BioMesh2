#pragma once

#include "biomesh2/BoundingBox.hpp"
#include "biomesh2/Atom.hpp"
#include <vector>
#include <memory>
#include <functional>
#include <array>
#include <atomic>
#include <mutex>
#include <unordered_set>
#include <unordered_map>

namespace biomesh2 {

// Forward declarations
class Octree;

/**
 * @brief Refinement criteria for adaptive octree subdivision
 */
enum class RefinementCriterion {
    NONE,                    // No specific criterion
    ATOM_DENSITY,           // Based on number of atoms in cell
    GEOMETRIC_FEATURE,      // Based on geometric features
    BOUNDARY_PROXIMITY,     // Based on proximity to molecular boundaries
    USER_DEFINED,           // User-defined criterion
    ERROR_INDICATOR         // Based on solution error indicators
};

/**
 * @brief Neighbor direction enumeration for efficient neighbor finding
 */
enum class NeighborDirection {
    // Face neighbors (6 directions)
    FACE_X_NEG = 0, FACE_X_POS = 1,
    FACE_Y_NEG = 2, FACE_Y_POS = 3,
    FACE_Z_NEG = 4, FACE_Z_POS = 5,
    
    // Edge neighbors (12 directions)
    EDGE_XY_NEG_NEG = 6, EDGE_XY_NEG_POS = 7,
    EDGE_XY_POS_NEG = 8, EDGE_XY_POS_POS = 9,
    EDGE_XZ_NEG_NEG = 10, EDGE_XZ_NEG_POS = 11,
    EDGE_XZ_POS_NEG = 12, EDGE_XZ_POS_POS = 13,
    EDGE_YZ_NEG_NEG = 14, EDGE_YZ_NEG_POS = 15,
    EDGE_YZ_POS_NEG = 16, EDGE_YZ_POS_POS = 17,
    
    // Vertex neighbors (8 directions)
    VERTEX_NEG_NEG_NEG = 18, VERTEX_NEG_NEG_POS = 19,
    VERTEX_NEG_POS_NEG = 20, VERTEX_NEG_POS_POS = 21,
    VERTEX_POS_NEG_NEG = 22, VERTEX_POS_NEG_POS = 23,
    VERTEX_POS_POS_NEG = 24, VERTEX_POS_POS_POS = 25
};

/**
 * @brief Enhanced octree node with advanced features for FEM-ready adaptive octree
 */
class OctreeNode {
public:
    // Basic geometric properties
    Point3D min;         // Minimum corner of the node's bounding box
    Point3D max;         // Maximum corner of the node's bounding box  
    Point3D center;      // Center point of the node
    Point3D halfSize;    // Half-size in each dimension
    
    // Tree structure
    int depth;           // Depth level in the tree (root = 0)
    int refinementLevel; // Refinement level (may differ from depth for adaptive refinement)
    bool isLeaf;         // True if this node has no children
    std::array<std::unique_ptr<OctreeNode>, 8> children; // 8 octant children
    
    // Enhanced tree navigation
    OctreeNode* parent;  // Parent node (nullptr for root)
    size_t nodeId;       // Unique node identifier
    int childIndex;      // Index in parent's children array (-1 for root)
    
    // Refinement and adaptation
    RefinementCriterion refinementCriterion; // Why this node was refined
    double refinementValue;                  // Quantitative refinement metric
    bool needsRefinement;                   // Flag for dynamic refinement
    bool needsCoarsening;                   // Flag for dynamic coarsening
    std::atomic<bool> isBeingProcessed;     // Thread safety flag
    
    // Content and analysis
    std::vector<size_t> containedAtomIds;   // IDs of atoms within this node
    double atomDensity;                     // Atoms per unit volume
    bool intersectsBoundary;                // True if node intersects molecular boundary
    double boundaryDistance;                // Distance to nearest molecular boundary
    
    // Neighbor management (for 2:1 balance constraint)
    std::unordered_map<NeighborDirection, std::vector<OctreeNode*>> neighbors;
    bool hasViolated2To1Constraint;         // Flag for constraint violation
    
    /**
     * @brief Constructor for octree node
     * @param minPoint Minimum corner coordinates
     * @param maxPoint Maximum corner coordinates
     * @param nodeDepth Depth level of this node
     * @param parentNode Parent node pointer (nullptr for root)
     * @param childIdx Index in parent's children array
     */
    OctreeNode(const Point3D& minPoint, const Point3D& maxPoint, 
               int nodeDepth = 0, OctreeNode* parentNode = nullptr, int childIdx = -1);
    
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
    
    /**
     * @brief Get the surface area of this node
     * @return Surface area in square units
     */
    double getSurfaceArea() const;
    
    /**
     * @brief Check if this node intersects with an atom
     * @param atom Atom to check intersection with
     * @return true if atom intersects with node bounds
     */
    bool intersectsAtom(const Atom& atom) const;
    
    /**
     * @brief Calculate atom density within this node
     * @param atoms Vector of all atoms
     * @return Number of atoms per unit volume
     */
    double calculateAtomDensity(const std::vector<std::unique_ptr<Atom>>& atoms);
    
    /**
     * @brief Update contained atoms list
     * @param atoms Vector of all atoms
     */
    void updateContainedAtoms(const std::vector<std::unique_ptr<Atom>>& atoms);
    
    /**
     * @brief Check if node needs refinement based on criteria
     * @param criterion Refinement criterion to check
     * @param threshold Threshold value for refinement
     * @param atoms Vector of atoms for density calculations
     * @return true if node should be refined
     */
    bool shouldRefine(RefinementCriterion criterion, double threshold,
                     const std::vector<std::unique_ptr<Atom>>& atoms = {}) const;
    
    /**
     * @brief Add neighbor in specified direction
     * @param direction Neighbor direction
     * @param neighbor Neighboring node
     */
    void addNeighbor(NeighborDirection direction, OctreeNode* neighbor);
    
    /**
     * @brief Get neighbors in specified direction
     * @param direction Neighbor direction
     * @return Vector of neighboring nodes
     */
    const std::vector<OctreeNode*>& getNeighbors(NeighborDirection direction) const;
    
    /**
     * @brief Check if this node violates 2:1 balance constraint
     * @return true if constraint is violated
     */
    bool violates2To1Constraint() const;
    
    /**
     * @brief Get minimum edge length of this node
     * @return Minimum edge length
     */
    double getMinEdgeLength() const;
    
    /**
     * @brief Get maximum edge length of this node
     * @return Maximum edge length
     */
    double getMaxEdgeLength() const;
    
    /**
     * @brief Thread-safe setter for processing flag
     * @param processing Processing state
     */
    void setProcessing(bool processing);
    
    /**
     * @brief Thread-safe getter for processing flag
     * @return Current processing state
     */
    bool getProcessing() const;

private:
    static size_t nextNodeId_;              // Static counter for unique node IDs
    static std::mutex nodeIdMutex_;         // Mutex for thread-safe ID generation
    
    /**
     * @brief Generate unique node ID
     * @return Unique node identifier
     */
    static size_t generateNodeId();
};

/**
 * @brief Advanced adaptive octree data structure for FEM-ready 3D space partitioning
 * 
 * Features:
 * - Adaptive refinement based on multiple criteria
 * - 2:1 balance constraint enforcement
 * - Efficient neighbor finding across refinement levels
 * - Boundary conformity for molecular structures
 * - Thread-safe operations for parallel processing
 * - Dynamic refinement and coarsening
 * - Memory-optimized lazy child allocation
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
     * @brief Basic subdivision with termination conditions
     * @param maxDepth Maximum depth to subdivide (default: 8)
     * @param minCellSize Minimum cell size before stopping subdivision (default: 0.001)
     * @param occupancyCheck Optional function to check if cell should be subdivided based on content
     */
    void subdivide(int maxDepth = 8, 
                   double minCellSize = 0.001,
                   std::function<bool(const OctreeNode&)> occupancyCheck = nullptr);
    
    /**
     * @brief Advanced adaptive subdivision with multiple criteria
     * @param atoms Vector of atoms for density and boundary calculations
     * @param maxDepth Maximum depth to subdivide
     * @param minCellSize Minimum cell size before stopping subdivision
     * @param criterion Primary refinement criterion
     * @param threshold Threshold value for refinement criterion
     * @param enforce2To1 Whether to enforce 2:1 balance constraint
     * @param enableBoundaryRefinement Whether to refine near molecular boundaries
     */
    void adaptiveSubdivide(const std::vector<std::unique_ptr<Atom>>& atoms,
                          int maxDepth = 8,
                          double minCellSize = 0.001,
                          RefinementCriterion criterion = RefinementCriterion::ATOM_DENSITY,
                          double threshold = 1.0,
                          bool enforce2To1 = true,
                          bool enableBoundaryRefinement = true);
    
    /**
     * @brief Enforce 2:1 balance constraint across the tree
     * @param atoms Vector of atoms (needed for refinement criteria)
     */
    void enforce2To1Balance(const std::vector<std::unique_ptr<Atom>>& atoms = {});
    
    /**
     * @brief Find all neighbors of a given node
     * @param node Node to find neighbors for
     * @param direction Specific direction (optional, finds all if not specified)
     * @return Vector of neighboring nodes
     */
    std::vector<OctreeNode*> findNeighbors(OctreeNode* node, 
                                          NeighborDirection direction = NeighborDirection::FACE_X_NEG) const;
    
    /**
     * @brief Find face neighbors (6 directions)
     * @param node Node to find neighbors for
     * @return Map of face neighbors by direction
     */
    std::unordered_map<NeighborDirection, std::vector<OctreeNode*>> 
    findFaceNeighbors(OctreeNode* node) const;
    
    /**
     * @brief Dynamically refine a specific node
     * @param node Node to refine
     * @param atoms Vector of atoms for refinement criteria
     * @param criterion Refinement criterion to use
     * @param threshold Threshold for refinement
     * @return true if refinement was successful
     */
    bool refineNode(OctreeNode* node, 
                   const std::vector<std::unique_ptr<Atom>>& atoms,
                   RefinementCriterion criterion = RefinementCriterion::ATOM_DENSITY,
                   double threshold = 1.0);
    
    /**
     * @brief Dynamically coarsen a node and its siblings
     * @param node Node to consider for coarsening
     * @param atoms Vector of atoms for coarsening criteria
     * @return true if coarsening was successful
     */
    bool coarsenNode(OctreeNode* node, 
                    const std::vector<std::unique_ptr<Atom>>& atoms);
    
    /**
     * @brief Update neighbor connectivity for entire tree
     */
    void updateNeighborConnectivity();
    
    /**
     * @brief Get all nodes at a specific refinement level
     * @param level Refinement level
     * @return Vector of nodes at the specified level
     */
    std::vector<OctreeNode*> getNodesAtLevel(int level) const;
    
    /**
     * @brief Get all leaf nodes suitable for FEM mesh generation
     * @return Vector of leaf nodes
     */
    std::vector<OctreeNode*> getFEMReadyLeaves() const;
    
    /**
     * @brief Validate tree structure and constraints
     * @return Vector of validation error messages (empty if valid)
     */
    std::vector<std::string> validateTree() const;
    
    /**
     * @brief Get the root node of the octree
     * @return Reference to root node
     */
    const OctreeNode& getRoot() const { return *root_; }
    
    /**
     * @brief Get mutable root node (for advanced operations)
     * @return Pointer to root node
     */
    OctreeNode* getRootPtr() { return root_.get(); }
    
    /**
     * @brief Print the tree structure in a readable format
     * @param showDepth Maximum depth to print (default: all)
     * @param showNeighbors Whether to show neighbor information
     */
    void printTree(int showDepth = -1, bool showNeighbors = false) const;
    
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
     * @brief Get number of nodes at each depth level
     * @return Vector where index is depth and value is node count
     */
    std::vector<size_t> getNodeCountByDepth() const;
    
    /**
     * @brief Get detailed tree statistics
     * @return Map of statistic name to value
     */
    std::unordered_map<std::string, double> getTreeStatistics() const;
    
    /**
     * @brief Find the leaf node that contains a given point
     * @param point Point to locate
     * @return Pointer to containing leaf node, nullptr if point is outside domain
     */
    const OctreeNode* findLeaf(const Point3D& point) const;
    
    /**
     * @brief Find node by unique ID
     * @param nodeId Unique node identifier
     * @return Pointer to node, nullptr if not found
     */
    OctreeNode* findNodeById(size_t nodeId) const;
    
    /**
     * @brief Check if tree satisfies 2:1 balance constraint
     * @return true if balanced, false otherwise
     */
    bool is2To1Balanced() const;
    
    /**
     * @brief Get memory usage statistics
     * @return Memory usage in bytes
     */
    size_t getMemoryUsage() const;

private:
    std::unique_ptr<OctreeNode> root_;
    mutable std::mutex treeMutex_;           // For thread-safe operations
    std::unordered_map<size_t, OctreeNode*> nodeRegistry_; // Fast node lookup by ID
    
    // Configuration
    bool enableThreadSafety_;
    bool enableMemoryOptimization_;
    double defaultRefinementThreshold_;
    
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
     * @brief Adaptive recursive subdivision with multiple criteria
     * @param node Node to subdivide
     * @param atoms Vector of atoms
     * @param maxDepth Maximum depth allowed
     * @param minCellSize Minimum cell size allowed
     * @param criterion Refinement criterion
     * @param threshold Refinement threshold
     */
    void adaptiveSubdivideRecursive(OctreeNode* node,
                                   const std::vector<std::unique_ptr<Atom>>& atoms,
                                   int maxDepth,
                                   double minCellSize,
                                   RefinementCriterion criterion,
                                   double threshold);
    
    /**
     * @brief Create child nodes with lazy allocation
     * @param parent Parent node
     * @param childIndex Index of child to create (or -1 for all)
     * @return Pointer to created child (or nullptr if childIndex == -1)
     */
    OctreeNode* createChild(OctreeNode* parent, int childIndex = -1);
    
    /**
     * @brief Create all 8 children for a parent node
     * @param parent Parent node
     */
    void createAllChildren(OctreeNode* parent);
    
    /**
     * @brief Update atom information for a node and its subtree
     * @param node Root node to update
     * @param atoms Vector of atoms
     */
    void updateNodeAtomInfo(OctreeNode* node, const std::vector<std::unique_ptr<Atom>>& atoms);
    
    /**
     * @brief Perform boundary refinement
     * @param atoms Vector of atoms
     * @param maxDepth Maximum depth allowed
     * @param minCellSize Minimum cell size
     */
    void performBoundaryRefinement(const std::vector<std::unique_ptr<Atom>>& atoms, 
                                  int maxDepth, double minCellSize);
    
    /**
     * @brief Collect nodes that violate 2:1 constraint
     * @param node Root node to search
     * @param violatingNodes Vector to store violating nodes
     */
    void collectViolatingNodes(OctreeNode* node, std::vector<OctreeNode*>& violatingNodes);
    
    /**
     * @brief Collect nodes that intersect boundaries
     * @param node Root node to search
     * @param boundaryNodes Vector to store boundary nodes
     */
    void collectBoundaryNodes(OctreeNode* node, std::vector<OctreeNode*>& boundaryNodes);
    
    /**
     * @brief Update neighbors after node subdivision
     * @param parent Parent node that was subdivided
     */
    void updateNeighborsAfterSubdivision(OctreeNode* parent);
    
    /**
     * @brief Recursive neighbor finding implementation
     * @param node Current node
     * @param target Target point or bounds
     * @param direction Search direction
     * @param level Search level
     * @return Vector of neighboring nodes
     */
    std::vector<OctreeNode*> findNeighborsRecursive(OctreeNode* node,
                                                   const Point3D& target,
                                                   NeighborDirection direction,
                                                   int level) const;
    
    /**
     * @brief Check and enforce 2:1 constraint for a specific node
     * @param node Node to check
     * @param atoms Vector of atoms (for refinement if needed)
     * @return true if constraint was enforced
     */
    bool constrain2To1ForNode(OctreeNode* node, 
                             const std::vector<std::unique_ptr<Atom>>& atoms);
    
    /**
     * @brief Clear neighbor connectivity for entire subtree
     * @param node Root node
     */
    void clearNeighborConnectivity(OctreeNode* node);
    
    /**
     * @brief Build neighbor connectivity for entire subtree
     * @param node Root node
     */
    void buildNeighborConnectivity(OctreeNode* node);
    
    /**
     * @brief Collect nodes at specific level
     * @param node Root node
     * @param level Target level
     * @param nodes Vector to store nodes
     */
    void collectNodesAtLevel(OctreeNode* node, int level, std::vector<OctreeNode*>& nodes) const;
    
    /**
     * @brief Collect FEM-ready leaves
     * @param node Root node
     * @param femLeaves Vector to store FEM-ready leaves
     */
    void collectFEMReadyLeaves(OctreeNode* node, std::vector<OctreeNode*>& femLeaves) const;
    
    /**
     * @brief Validate tree structure recursively
     * @param node Root node
     * @param errors Vector to store error messages
     */
    void validateTreeStructure(OctreeNode* node, std::vector<std::string>& errors) const;
    
    /**
     * @brief Validate neighbor connectivity
     * @param node Root node
     * @param errors Vector to store error messages
     */
    void validateNeighborConnectivity(OctreeNode* node, std::vector<std::string>& errors) const;
    
    /**
     * @brief Check balance constraint recursively
     * @param node Root node
     * @return true if balanced
     */
    bool checkBalanceRecursive(OctreeNode* node) const;
    
    /**
     * @brief Calculate memory usage recursively
     * @param node Root node
     * @return Memory usage in bytes
     */
    size_t calculateMemoryUsage(OctreeNode* node) const;
    
    /**
     * @brief Get maximum depth of tree
     * @param node Root node
     * @return Maximum depth
     */
    int getMaxDepth(OctreeNode* node) const;
    
    /**
     * @brief Get minimum leaf size
     * @param node Root node
     * @return Minimum leaf size
     */
    double getMinLeafSize(OctreeNode* node) const;
    
    /**
     * @brief Get maximum leaf size
     * @param node Root node
     * @return Maximum leaf size
     */
    double getMaxLeafSize(OctreeNode* node) const;
    
    /**
     * @brief Get average atom density
     * @param node Root node
     * @return Average atom density
     */
    double getAverageAtomDensity(OctreeNode* node) const;
    
    /**
     * @brief Calculate average atom density recursively
     * @param node Current node
     * @param totalDensity Accumulated density
     * @param leafCount Number of leaves processed
     */
    void calculateAverageAtomDensityRecursive(OctreeNode* node, double& totalDensity, size_t& leafCount) const;
    
    /**
     * @brief Print tree structure recursively with enhanced information
     * @param node Current node to print
     * @param prefix String prefix for indentation
     * @param maxDepth Maximum depth to print
     * @param showNeighbors Whether to show neighbor information
     */
    void printTreeRecursive(const OctreeNode* node, const std::string& prefix, 
                           int maxDepth, bool showNeighbors) const;
    
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
     * @brief Count nodes by depth recursively
     * @param node Current node
     * @param counts Vector to store counts by depth
     */
    void countNodesByDepth(const OctreeNode* node, std::vector<size_t>& counts) const;
    
    /**
     * @brief Find leaf containing point recursively
     * @param node Current node to search
     * @param point Point to locate
     * @return Pointer to containing leaf node or nullptr
     */
    const OctreeNode* findLeafRecursive(const OctreeNode* node, const Point3D& point) const;
    
    /**
     * @brief Find node by ID recursively
     * @param node Current node to search
     * @param nodeId Target node ID
     * @return Pointer to node or nullptr
     */
    OctreeNode* findNodeByIdRecursive(OctreeNode* node, size_t nodeId) const;
    
    /**
     * @brief Register node in the node registry
     * @param node Node to register
     */
    void registerNode(OctreeNode* node);
    
    /**
     * @brief Unregister node from the node registry
     * @param node Node to unregister
     */
    void unregisterNode(OctreeNode* node);
    
    /**
     * @brief Thread-safe tree operation wrapper
     * @param operation Function to execute with tree lock
     */
    template<typename Func>
    auto lockAndExecute(Func&& operation) const -> decltype(operation()) {
        if (enableThreadSafety_) {
            std::lock_guard<std::mutex> lock(treeMutex_);
            return operation();
        } else {
            return operation();
        }
    }
};

} // namespace biomesh2