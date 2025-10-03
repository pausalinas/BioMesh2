#pragma once

#include "biomesh2/VoxelGrid.hpp"
#include "biomesh2/HexMesh.hpp"
#include <vector>
#include <array>
#include <unordered_map>

namespace biomesh2 {

/**
 * @brief Voxel grid to hexahedral mesh generator
 * 
 * Generates a hexahedral finite element mesh from a uniform voxel grid.
 * Only occupied voxels (those intersecting with atom spheres) are converted
 * to hexahedral elements.
 */
class VoxelMeshGenerator {
public:
    /**
     * @brief Generate hexahedral mesh from voxel grid occupied voxels
     * @param voxelGrid The voxel grid to generate mesh from
     * @return HexMesh containing nodes and element connectivity
     */
    static HexMesh generateHexMesh(const VoxelGrid& voxelGrid);

private:
    /**
     * @brief Compute the 8 corner nodes for a given voxel
     * @param voxel The voxel to compute corners for
     * @return Array of 8 corner points in standard hexahedral ordering
     */
    static std::array<Point3D, 8> computeCornerNodes(const Voxel& voxel);
    
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
