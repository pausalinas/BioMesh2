#include "biomesh2/HexMesh.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

namespace biomesh2 {

bool HexMesh::exportToMsh(const std::string& filename) const {
    std::ofstream outFile(filename);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return false;
    }
    
    // Write Gmsh format header (version 2.2, ASCII)
    outFile << "$MeshFormat\n";
    outFile << "2.2 0 8\n";
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
        // Element type 5 = 8-node hexahedron
        // Format: elm-number elm-type number-of-tags < tags > node-number-list
        outFile << (i + 1) << " 5 0";  // 5 = hexahedron, 0 tags
        
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

std::string HexMesh::generateMeshFilename(const std::string& pdbFilePath, double resolution) {
    // Extract PDB code from filename
    // Example: "path/to/1ABC.pdb" -> "1ABC"
    // or "test_peptide.pdb" -> "test_peptide"
    
    size_t lastSlash = pdbFilePath.find_last_of("/\\");
    std::string filename = (lastSlash != std::string::npos) 
                          ? pdbFilePath.substr(lastSlash + 1) 
                          : pdbFilePath;
    
    // Remove .pdb extension if present
    size_t dotPos = filename.find_last_of(".");
    if (dotPos != std::string::npos) {
        filename = filename.substr(0, dotPos);
    }
    
    // Format resolution with one decimal place
    std::ostringstream resStream;
    resStream << std::fixed << std::setprecision(1) << resolution;
    std::string resStr = resStream.str();
    
    // Replace decimal point with underscore for filename compatibility
    std::replace(resStr.begin(), resStr.end(), '.', '_');
    
    // Generate filename: {PDB_code}{resolution}.msh
    return filename + resStr + ".msh";
}

std::string HexMesh::generateLogFilename(const std::string& pdbFilePath) {
    // Extract PDB code from filename
    size_t lastSlash = pdbFilePath.find_last_of("/\\");
    std::string filename = (lastSlash != std::string::npos) 
                          ? pdbFilePath.substr(lastSlash + 1) 
                          : pdbFilePath;
    
    // Remove .pdb extension if present
    size_t dotPos = filename.find_last_of(".");
    if (dotPos != std::string::npos) {
        filename = filename.substr(0, dotPos);
    }
    
    // Generate filename: {PDB_code}_output.log
    return filename + "_output.log";
}

} // namespace biomesh2
