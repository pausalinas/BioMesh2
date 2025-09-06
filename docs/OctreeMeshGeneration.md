# Octree to Hexahedral Mesh Generation

This module provides functionality to generate hexahedral finite element meshes from octree data structures. The implementation processes only the leaf nodes of the octree, creating one hexahedral element for each leaf cell.

## Features

- **Leaf Node Processing**: Traverses the octree and collects all leaf cells (nodes with no children)
- **Node Deduplication**: Automatically detects and reuses shared corner nodes between adjacent elements
- **Standard Element Ordering**: Uses standard hexahedral element node ordering compatible with finite element software
- **Memory Efficient**: Uses hash maps for O(1) node lookup and deduplication
- **GiD Export**: Export generated meshes to GiD format for use with finite element analysis software

## Usage

```cpp
#include "biomesh2/OctreeMeshGenerator.hpp"

// Create and subdivide an octree
biomesh2::Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
octree.subdivide(3); // Subdivide to depth 3

// Generate hexahedral mesh
biomesh2::HexMesh mesh = biomesh2::OctreeMeshGenerator::generateHexMesh(octree);

// Access mesh data
std::cout << "Nodes: " << mesh.getNodeCount() << "\n";
std::cout << "Elements: " << mesh.getElementCount() << "\n";

// Access individual nodes and elements
const biomesh2::Point3D& node0 = mesh.nodes[0];
const std::array<int, 8>& element0 = mesh.elements[0];

// Export mesh to GiD format
biomesh2::OctreeMeshGenerator::exportToGiD(mesh, "output_mesh.msh");
```

### GiD Export Usage

```cpp
#include "biomesh2/OctreeMeshGenerator.hpp"

// Generate mesh from octree
biomesh2::Octree octree(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
octree.subdivide(2);
biomesh2::HexMesh mesh = biomesh2::OctreeMeshGenerator::generateHexMesh(octree);

// Export to GiD format
try {
    biomesh2::OctreeMeshGenerator::exportToGiD(mesh, "my_mesh.msh");
    std::cout << "Successfully exported mesh to my_mesh.msh" << std::endl;
} catch (const std::exception& e) {
    std::cerr << "Export failed: " << e.what() << std::endl;
}
```

## Data Structures

### HexMesh
Contains the generated mesh data:
- `nodes`: Vector of unique 3D node coordinates
- `elements`: Vector of element connectivity arrays (8 node indices per hexahedron)

### Node Ordering
Each hexahedral element uses the following node ordering:
```
       7--------6
      /|       /|
     / |      / |
    4--------5  |
    |  3-----|--2
    | /      | /
    |/       |/
    0--------1
```

- Nodes 0-3: Bottom face (z = min)
- Nodes 4-7: Top face (z = max)
- Node ordering follows right-hand rule

## GiD Mesh Format

The exported GiD mesh files use the following format:

```
MESH dimension 3 Elemtype Hexahedra Nnode 8
coordinates
1 x1 y1 z1
2 x2 y2 z2
...
end coordinates
elements
1 1 2 3 4 5 6 7 8
2 9 10 11 12 13 14 15 16
...
end elements
```

Key features of the GiD export:
- Uses 1-based indexing for both nodes and elements
- Coordinates are written with 6 decimal places precision
- Element connectivity follows the standard hexahedral node ordering
- Compatible with GiD pre/post-processor and other FEM software

## Implementation Details

1. **Collect Leaf Nodes**: Recursively traverses the octree to find all leaf cells
2. **Compute Corner Nodes**: For each leaf cell, computes its 8 corner coordinates
3. **Assign Unique Indices**: Uses hash map to deduplicate shared nodes between adjacent elements
4. **Build Connectivity**: Creates element connectivity arrays mapping to unique node indices

The node deduplication uses a tolerance of 1e-12 for floating point coordinate comparison to handle numerical precision issues.