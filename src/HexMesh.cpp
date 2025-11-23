#include "biomesh2/HexMesh.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace biomesh2 {

// Gmsh format constants
namespace {
    constexpr const char* GMSH_VERSION = "2.2";
    constexpr int GMSH_FILE_TYPE = 0;  // ASCII format
    constexpr int GMSH_DATA_SIZE = 8;  // Size of double in bytes
    constexpr int GMSH_HEXAHEDRON_TYPE = 5;  // Element type for 8-node hexahedron
    constexpr int GMSH_NO_TAGS = 0;  // No tags for elements
}

bool HexMesh::exportToMsh(const std::string& filename) const {
    std::ofstream outFile(filename);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return false;
    }
    
    // Write Gmsh format header
    outFile << "$MeshFormat\n";
    outFile << GMSH_VERSION << " " << GMSH_FILE_TYPE << " " << GMSH_DATA_SIZE << "\n";
    outFile << "$EndMeshFormat\n";
    
    // Write nodes section
    outFile << "$Nodes\n";
    outFile << nodes.size() << "\n";
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        // Node numbering starts from 1 in Gmsh format
        outFile << (i + 1) << " " 
                << std::fixed << std::setprecision(10) 
                << nodes[i].x << " " 
                << nodes[i].y << " " 
                << nodes[i].z << "\n";
    }
    
    outFile << "$EndNodes\n";
    
    // Write elements section
    outFile << "$Elements\n";
    outFile << elements.size() << "\n";
    
    for (size_t i = 0; i < elements.size(); ++i) {
        // Element numbering starts from 1
        // Format: elm-number elm-type number-of-tags < tags > node-number-list
        outFile << (i + 1) << " " << GMSH_HEXAHEDRON_TYPE << " " << GMSH_NO_TAGS;
        
        // Write the 8 node indices (convert from 0-based to 1-based indexing)
        for (int j = 0; j < 8; ++j) {
            outFile << " " << (elements[i][j] + 1);
        }
        outFile << "\n";
    }
    
    outFile << "$EndElements\n";
    
    outFile.close();
    
    std::cout << "Successfully exported mesh to: " << filename << std::endl;
    std::cout << "  Nodes: " << nodes.size() << std::endl;
    std::cout << "  Elements: " << elements.size() << std::endl;
    
    return true;
}

std::string HexMesh::extractBaseName(const std::string& pdbFilePath) {
    // Extract filename from path
    size_t lastSlash = pdbFilePath.find_last_of("/\\");
    std::string filename = (lastSlash != std::string::npos) 
                          ? pdbFilePath.substr(lastSlash + 1) 
                          : pdbFilePath;
    
    // Remove .pdb extension if present (case-insensitive)
    const std::string pdbExt = ".pdb";
    if (filename.size() >= pdbExt.size()) {
        std::string ending = filename.substr(filename.size() - pdbExt.size());
        // Convert to lowercase for case-insensitive comparison
        std::transform(ending.begin(), ending.end(), ending.begin(), ::tolower);
        if (ending == pdbExt) {
            filename = filename.substr(0, filename.size() - pdbExt.size());
        }
    }
    
    return filename;
}

std::string HexMesh::generateMeshFilename(const std::string& pdbFilePath, double resolution) {
    // Extract base name from PDB file path
    std::string baseName = extractBaseName(pdbFilePath);
    
    // Format resolution with one decimal place
    std::ostringstream resStream;
    resStream << std::fixed << std::setprecision(1) << resolution;
    std::string resStr = resStream.str();
    
    // Replace decimal point with underscore for filename compatibility
    std::replace(resStr.begin(), resStr.end(), '.', '_');
    
    // Generate filename: {PDB_code}{resolution}.msh
    return baseName + resStr + ".msh";
}

std::string HexMesh::generateLogFilename(const std::string& pdbFilePath) {
    // Extract base name from PDB file path
    std::string baseName = extractBaseName(pdbFilePath);
    
    // Generate filename: {PDB_code}_output.log
    return baseName + "_output.log";
}

} // namespace biomesh2
