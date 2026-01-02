#pragma once

#include "biomesh/HexMesh.hpp"
#include <string>

namespace biomesh {

/**
 * @brief VTK mesh file exporter
 *
 * Exports hexahedral meshes to the VTK legacy ASCII format (.vtk files)
 * using the UNSTRUCTURED_GRID dataset with cell type 12 (VTK_HEXAHEDRON).
 */
class VTKExporter {
public:
    /**
     * @brief Export hexahedral mesh to VTK format (.vtk file)
     * @param mesh The hexahedral mesh to export
     * @param filename Output file path (should end with .vtk)
     * @return true if export was successful, false otherwise
     */
    static bool exportToVTK(const HexMesh& mesh, const std::string& filename);
};

} // namespace biomesh
