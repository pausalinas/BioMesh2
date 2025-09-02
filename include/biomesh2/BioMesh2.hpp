#pragma once

/**
 * @file BioMesh2.hpp
 * @brief Main header file for the BioMesh2 C++ module
 * 
 * This module provides functionality for:
 * - Parsing PDB structure files
 * - Extracting and enriching atom information with physical properties
 * - Computing molecular bounding boxes
 * - Octree spatial partitioning and mesh generation
 * - Finite element stress analysis and von Mises stress calculations
 * - Statistical analysis for mutation scenario comparisons
 * 
 * @author BioMesh2 Team
 * @version 1.0.0
 */

#include "biomesh2/Atom.hpp"
#include "biomesh2/AtomicSpec.hpp"
#include "biomesh2/PDBParser.hpp"
#include "biomesh2/AtomBuilder.hpp"
#include "biomesh2/BoundingBox.hpp"
#include "biomesh2/Octree.hpp"
#include "biomesh2/OctreeMeshGenerator.hpp"
#include "biomesh2/StressAnalysis.hpp"
#include "biomesh2/StatisticalAnalysis.hpp"

namespace biomesh2 {

/**
 * @brief Main workflow function for processing PDB files
 * 
 * This convenience function combines all the steps:
 * 1. Parse PDB file
 * 2. Enrich atoms with physical properties
 * 3. Calculate bounding box
 * 
 * @param pdbFilename Path to PDB file
 * @param padding Additional padding for bounding box (default: 0.0)
 * @return Pair of (enriched atoms, bounding box)
 * @throws std::runtime_error if file cannot be processed
 */
std::pair<std::vector<std::unique_ptr<Atom>>, BoundingBox> 
processPDBFile(const std::string& pdbFilename, double padding = 0.0);

} // namespace biomesh2