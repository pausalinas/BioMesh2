#include "biomesh2/BioMesh2.hpp"
#include "biomesh2/OctreeMeshGenerator.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace biomesh2;

int main() {
    std::cout << "BioMesh2 Complete Pipeline: PDB to GiD Mesh Export\n";
    std::cout << "==================================================\n\n";
    
    try {
        // Step 1: Parse PDB data (using built-in test data)
        std::cout << "Step 1: Parsing molecular structure...\n";
        std::string testPdbContent = 
            "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N  \n"
            "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00           C  \n"
            "ATOM      3  C   ALA A   1      17.900  17.100   8.900  1.00 20.00           C  \n"
            "ATOM      4  O   ALA A   1      18.100  18.200   8.500  1.00 20.00           O  \n"
            "ATOM      5  CB  ALA A   1      17.900  17.100   8.900  1.00 20.00           C  \n"
            "END                                                                             \n";
        
        auto basicAtoms = PDBParser::parsePDBContent(testPdbContent);
        std::cout << "   Parsed " << basicAtoms.size() << " atoms from PDB data\n";
        
        // Step 2: Enrich atoms with properties
        std::cout << "\nStep 2: Enriching atoms with physical properties...\n";
        AtomBuilder builder;
        auto enrichedAtoms = builder.buildAtoms(basicAtoms);
        std::cout << "   Enriched " << enrichedAtoms.size() << " atoms with radii and masses\n";
        
        // Step 3: Calculate bounding box
        std::cout << "\nStep 3: Computing molecular bounding box...\n";
        double padding = 2.0; // 2 Angstrom padding
        BoundingBox bbox(enrichedAtoms, padding);
        
        std::cout << "   Bounding box with " << padding << "Å padding:\n";
        std::cout << "   Min: (" << std::fixed << std::setprecision(3) 
                  << bbox.getMin().x << ", " << bbox.getMin().y << ", " << bbox.getMin().z << ")\n";
        std::cout << "   Max: (" << bbox.getMax().x << ", " << bbox.getMax().y << ", " << bbox.getMax().z << ")\n";
        std::cout << "   Volume: " << bbox.getVolume() << " ų\n";
        
        // Step 4: Create octree from bounding box
        std::cout << "\nStep 4: Creating octree spatial structure...\n";
        Point3D min = bbox.getMin();
        Point3D max = bbox.getMax();
        Octree octree(min.x, min.y, min.z, 
                     max.x - min.x, max.y - min.y, max.z - min.z);
        
        // Subdivide octree to create a reasonable mesh resolution
        octree.subdivide(2); // Creates 64 leaf cells
        std::cout << "   Created octree with " << octree.getNodeCount() << " nodes\n";
        std::cout << "   Number of leaf cells: " << octree.getLeafCount() << "\n";
        
        // Step 5: Generate hexahedral mesh
        std::cout << "\nStep 5: Generating hexahedral mesh...\n";
        HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
        std::cout << "   Generated mesh with:\n";
        std::cout << "   - Nodes: " << mesh.getNodeCount() << "\n";
        std::cout << "   - Elements: " << mesh.getElementCount() << "\n";
        
        // Step 6: Export to GiD format
        std::cout << "\nStep 6: Exporting mesh to GiD format...\n";
        std::string meshFilename = "biomolecular_mesh.msh";
        OctreeMeshGenerator::exportToGiD(mesh, meshFilename);
        std::cout << "   ✓ Successfully exported mesh to " << meshFilename << "\n";
        
        // Step 7: Show sample of the exported file
        std::cout << "\nStep 7: Sample of exported GiD mesh file:\n";
        std::ifstream file(meshFilename);
        if (file.is_open()) {
            std::string line;
            int lineCount = 0;
            while (std::getline(file, line) && lineCount < 20) {
                std::cout << "   " << line << "\n";
                lineCount++;
            }
            if (lineCount >= 20) {
                std::cout << "   ... (file continues with remaining " 
                          << (mesh.getNodeCount() + mesh.getElementCount()) - 10 
                          << " data lines)\n";
            }
            file.close();
        }
        
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "SUCCESS: Complete pipeline executed successfully!\n";
        std::cout << "\nWorkflow Summary:\n";
        std::cout << "  PDB atoms parsed    → " << basicAtoms.size() << " atoms\n";
        std::cout << "  Properties enriched → " << enrichedAtoms.size() << " atoms with radii/masses\n";
        std::cout << "  Bounding box        → " << std::setprecision(1) << bbox.getVolume() << " ų volume\n";
        std::cout << "  Octree subdivision  → " << octree.getLeafCount() << " leaf cells\n";
        std::cout << "  Mesh generation     → " << mesh.getNodeCount() << " nodes, " << mesh.getElementCount() << " elements\n";
        std::cout << "  GiD export          → " << meshFilename << " (ready for FEM analysis)\n";
        std::cout << "\nThe exported mesh can now be imported into GiD, Abaqus, ANSYS, or other FEM software.\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error in pipeline: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}