#include "biomesh/EmptyVoxelMeshGenerator.hpp"
#include "biomesh/GiDExporter.hpp"
#include "biomesh/VTKExporter.hpp"
#include "biomesh/VoxelGrid.hpp"
#include "biomesh/PDBParser.hpp"
#include "biomesh/AtomBuilder.hpp"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace biomesh;

void printUsage(const char* programName) {
    std::cout << "\nUsage: " << programName << " <pdb_file> <voxel_size> <output_file> [padding]\n\n";
    std::cout << "Arguments:\n";
    std::cout << "  pdb_file     : Path to input PDB file\n";
    std::cout << "  voxel_size   : Edge length of voxels in Angstroms (e.g., 1.0)\n";
    std::cout << "  output_file  : Output mesh path; use .vtk for VTK or .msh/.gid for GiD\n";
    std::cout << "  padding      : Optional padding around bounding box in Angstroms (default: 2.0)\n\n";
    std::cout << "Example:\n";
    std::cout << "  " << programName << " protein.pdb 1.0 empty_mesh.vtk 2.0\n\n";
    std::cout << "Note: Empty voxel meshes can be very large (typically 95-99% of total voxels).\n";
    std::cout << "      Use larger voxel sizes for initial testing.\n\n";
}

std::string toLower(const std::string& value) {
    std::string lowered = value;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(), [](unsigned char c) { return std::tolower(c); });
    return lowered;
}

std::string detectFormat(const std::filesystem::path& outputPath) {
    const std::string ext = toLower(outputPath.extension().string());
    if (ext == ".vtk") {
        return "vtk";
    }
    if (ext == ".msh" || ext == ".gid") {
        return "gid";
    }
    return "gid";
}

std::string formatResolution(double value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << value;
    std::string res = oss.str();
    res.erase(res.find_last_not_of('0') + 1);
    if (!res.empty() && res.back() == '.') {
        res.pop_back();
    }
    return res;
}

std::string pickExtension(const std::filesystem::path& outputPath, const std::string& format) {
    std::string ext = toLower(outputPath.extension().string());
    if (!ext.empty() && ext[0] == '.') {
        ext.erase(0, 1);
    }
    if (!ext.empty()) {
        return ext;
    }
    return format;
}

std::string buildOutputFilename(const std::string& requestedOutput,
                                const std::string& pdbPath,
                                double voxelSize,
                                const std::string& format,
                                const std::string& extension) {
    namespace fs = std::filesystem;
    fs::path basePath(requestedOutput);
    std::string pdbId = fs::path(pdbPath).stem().string();
    if (pdbId.empty()) {
        pdbId = "pdb";
    }

    std::string resolutionTag = formatResolution(voxelSize);
    std::string stem = basePath.stem().string();
    if (stem.empty()) {
        stem = "mesh";
    }

    std::ostringstream filename;
    filename << stem << "_" << pdbId << "_res" << resolutionTag << "_" << format;

    const std::string& finalExt = extension.empty() ? format : extension;
    fs::path finalPath = basePath.parent_path() / (filename.str() + "." + finalExt);
    return finalPath.string();
}

int main(int argc, char* argv[]) {
    std::cout << "BioMesh - Empty Voxel Mesh Generator with GiD/VTK Export\n";
    std::cout << "======================================================\n\n";
    
    // Check arguments
    if (argc < 4) {
        std::cerr << "Error: Insufficient arguments\n";
        printUsage(argv[0]);
        return 1;
    }
    
    std::string pdbFile = argv[1];
    double voxelSize = 0.0;
    std::string outputFile = argv[3];
    double padding = 2.0; // Default padding
    
    // Parse voxel size
    try {
        voxelSize = std::stod(argv[2]);
        if (voxelSize <= 0.0) {
            std::cerr << "Error: Voxel size must be positive\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid voxel size: " << argv[2] << "\n";
        return 1;
    }
    
    // Parse optional padding
    if (argc > 4) {
        try {
            padding = std::stod(argv[4]);
            if (padding < 0.0) {
                std::cerr << "Error: Padding cannot be negative\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid padding value: " << argv[4] << "\n";
            return 1;
        }
    }
    
    try {
        // Step 1: Parse PDB file
        std::cout << "Loading PDB file: " << pdbFile << "\n";
        auto basicAtoms = PDBParser::parsePDBFile(pdbFile);
        std::cout << "  Loaded " << basicAtoms.size() << " atoms\n\n";
        
        // Step 2: Enrich atoms with physical properties
        std::cout << "Enriching atoms with physical properties...\n";
        AtomBuilder builder;
        auto enrichedAtoms = builder.buildAtoms(basicAtoms);
        std::cout << "  Enriched " << enrichedAtoms.size() << " atoms\n\n";
        
        // Step 3: Create voxel grid
        std::cout << "Creating voxel grid...\n";
        std::cout << "  Voxel size: " << voxelSize << " Å\n";
        std::cout << "  Padding: " << padding << " Å\n";
        VoxelGrid voxelGrid(enrichedAtoms, voxelSize, padding);
        
        std::cout << "\n";
        voxelGrid.printStatistics();
        std::cout << "\n";
        
        // Step 4: Generate hexahedral mesh from empty voxels
        std::cout << "Generating hexahedral mesh from empty voxels...\n";
        HexMesh mesh = EmptyVoxelMeshGenerator::generateHexMesh(voxelGrid);
        
        std::cout << "  Generated mesh:\n";
        std::cout << "    Nodes: " << mesh.getNodeCount() << "\n";
        std::cout << "    Elements: " << mesh.getElementCount() << "\n";
        
        // Calculate mesh statistics
        int emptyVoxelCount = voxelGrid.getEmptyVoxelCount();
        int theoreticalNodes = emptyVoxelCount * 8;
        double efficiency = 0.0;
        if (theoreticalNodes > 0) {
            efficiency = (1.0 - (double)mesh.getNodeCount() / theoreticalNodes) * 100.0;
        }
        std::cout << "    Node sharing efficiency: " << std::fixed << std::setprecision(1) 
                  << efficiency << "%\n\n";
        
        // Warn about large meshes
        if (mesh.getElementCount() > 100000) {
            std::cout << "WARNING: Large mesh detected (" << mesh.getElementCount() 
                      << " elements). File may be large.\n\n";
        }
        
        // Step 5: Export to selected format with metadata in the filename
        std::filesystem::path requestedPath(outputFile);
        std::string format = detectFormat(requestedPath);
        std::string extension = pickExtension(requestedPath, format);
        std::string finalOutput = buildOutputFilename(outputFile, pdbFile, voxelSize, format, extension);

        std::cout << "Exporting to " << format << " format: " << finalOutput << "\n";
        bool success = false;
        if (format == "vtk") {
            success = VTKExporter::exportToVTK(mesh, finalOutput);
        } else {
            success = GiDExporter::exportToGiD(mesh, finalOutput);
        }
        
        if (success) {
            std::cout << "  Export successful!\n";
            std::cout << "\nMesh file written to: " << finalOutput << "\n";
            std::cout << "You can now open this file in GiD/VTK-compatible FEM/CFD software.\n";
            return 0;
        } else {
            std::cerr << "  Export failed!\n";
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
