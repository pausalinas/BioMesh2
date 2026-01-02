#include "biomesh/VTKExporter.hpp"
#include <fstream>
#include <iostream>

namespace biomesh {

bool VTKExporter::exportToVTK(const HexMesh& mesh, const std::string& filename) {
    if (mesh.nodes.empty() || mesh.elements.empty()) {
        std::cerr << "Error: Cannot export empty mesh to VTK format\n";
        return false;
    }

    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << "\n";
        return false;
    }

    outFile << "# vtk DataFile Version 3.0\n";
    outFile << "BioMesh HexMesh\n";
    outFile << "ASCII\n";
    outFile << "DATASET UNSTRUCTURED_GRID\n";

    // Points section
    outFile << "POINTS " << mesh.nodes.size() << " double\n";
    for (const auto& node : mesh.nodes) {
        outFile << node.x << " " << node.y << " " << node.z << "\n";
    }

    // Cells section
    const size_t cellCount = mesh.elements.size();
    const size_t indicesPerCell = 9; // 8 nodes + leading count
    outFile << "CELLS " << cellCount << " " << (cellCount * indicesPerCell) << "\n";
    for (const auto& element : mesh.elements) {
        outFile << "8";
        for (int idx : element) {
            outFile << " " << idx;
        }
        outFile << "\n";
    }

    // Cell types section (12 = VTK_HEXAHEDRON)
    outFile << "CELL_TYPES " << cellCount << "\n";
    for (size_t i = 0; i < cellCount; ++i) {
        outFile << "12\n";
    }

    outFile.close();
    return true;
}

} // namespace biomesh
