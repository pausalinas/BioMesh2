# BioMesh2 - C++ PDB Parser and Molecular Bounding Box Calculator

A modern C++ module for parsing PDB (Protein Data Bank) structure files, extracting atom information, enriching it with physical properties, computing molecular bounding boxes, and generating hexahedral meshes using both octree and voxelization approaches.

## Features

- **PDB File Parsing**: Read PDB files and extract atom coordinates and element information
- **Atomic Property Enrichment**: Automatically assign atomic radii and masses based on element type
- **Bounding Box Calculation**: Compute 3D bounding boxes that encompass all atoms including their van der Waals radii
- **Octree Spatial Partitioning**: Efficient adaptive 3D space subdivision for spatial queries and geometric operations
- **Uniform Voxelization**: Tessellate molecular domains into regular cubic voxels for mesh generation
- **Hexahedral Mesh Generation**: Generate finite element meshes using both octree and voxelization approaches
- **Modern C++ Design**: Uses RAII, STL containers, smart pointers, and follows best practices
- **Unit Tested**: Comprehensive GoogleTest test suite with 57 tests
- **Extensible**: Easy to add new elements and their properties

## Project Structure

```
BioMesh2/
├── include/biomesh2/          # Header files
│   ├── Atom.hpp              # Enhanced Atom class with physical properties
│   ├── AtomicSpec.hpp        # Atomic specifications database
│   ├── PDBParser.hpp         # PDB file parsing functionality
│   ├── AtomBuilder.hpp       # Atom property enrichment
│   ├── BoundingBox.hpp       # Molecular bounding box calculation
│   ├── Octree.hpp            # Octree spatial partitioning data structure
│   ├── OctreeMeshGenerator.hpp  # Octree to hexahedral mesh converter
│   ├── VoxelGrid.hpp         # Uniform voxelization data structure
│   ├── VoxelMeshGenerator.hpp   # Voxel grid to hexahedral mesh converter
│   └── BioMesh2.hpp          # Main header with convenience functions
├── src/                      # Source files
├── tests/                    # GoogleTest unit tests
├── examples/                 # Example usage
│   ├── main.cpp             # Basic PDB parsing and bounding box demo
│   ├── octree_demo.cpp      # Octree functionality demonstration
│   └── voxel_demo.cpp       # Voxelization functionality demonstration
├── data/                     # Sample PDB files
├── docs/                     # Documentation
│   ├── OctreeMeshGeneration.md      # Octree meshing documentation
│   └── VoxelizationMeshGeneration.md # Voxelization meshing documentation
└── CMakeLists.txt           # CMake build configuration
```

## Building

### Prerequisites

- C++17 compatible compiler (GCC, Clang, or MSVC)
- CMake 3.14 or later
- GoogleTest (optional, for running tests)

### Build Instructions

```bash
# Clone the repository
git clone https://github.com/pausalinas/BioMesh2.git
cd BioMesh2

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build the project
make -j4

# Run tests (if GoogleTest is available)
make test

# Run example
./biomesh2_example

# Run octree demonstration
./octree_demo

# Run voxelization demonstration
./voxel_demo
```

## Usage

### Basic Usage

```cpp
#include "biomesh2/BioMesh2.hpp"
using namespace biomesh2;

// Process a PDB file in one step
auto [atoms, boundingBox] = processPDBFile("protein.pdb", 2.0); // 2.0 Å padding

// Access bounding box information
std::cout << "Volume: " << boundingBox.getVolume() << " Ų" << std::endl;
std::cout << "Dimensions: " << boundingBox.getDimensions().x << " × " 
          << boundingBox.getDimensions().y << " × " 
          << boundingBox.getDimensions().z << " Å" << std::endl;
```

### Step-by-Step Usage

```cpp
#include "biomesh2/BioMesh2.hpp"
using namespace biomesh2;

// Step 1: Parse PDB file
auto basicAtoms = PDBParser::parsePDBFile("protein.pdb");

// Step 2: Enrich with physical properties
AtomBuilder builder;
auto enrichedAtoms = builder.buildAtoms(basicAtoms);

// Step 3: Calculate bounding box
BoundingBox boundingBox(enrichedAtoms, 1.5); // 1.5 Å padding

// Access individual atom properties
for (const auto& atom : enrichedAtoms) {
    std::cout << "Element: " << atom->getChemicalElement()
              << ", Radius: " << atom->getAtomicRadius() << " Å"
              << ", Mass: " << atom->getAtomicMass() << " Da" << std::endl;
}
```

## Supported Elements

The module includes built-in atomic specifications for common elements found in biological molecules:

- **Main elements**: H, C, N, O, P, S
- **Halogens**: F, Cl, Br, I  
- **Metals**: Na, Mg, K, Ca, Fe, Zn
- **Others**: Se

Each element includes van der Waals radius and atomic mass. Additional elements can be easily added to the database.

## API Reference

### Atom Class

Enhanced atom structure with multiple constructors:

```cpp
// Constructor with element only
Atom atom("C");

// Constructor with element and radius
Atom atom("N", 1.55);

// Constructor with element, radius, and mass
Atom atom("O", 1.52, 15.999);
```

### AtomicSpecDatabase

Singleton database for atomic properties:

```cpp
auto& db = AtomicSpecDatabase::getInstance();
if (db.hasElement("C")) {
    const auto& spec = db.getSpec("C");
    std::cout << "Carbon radius: " << spec.radius << " Å" << std::endl;
}
```

### BoundingBox

Comprehensive bounding box functionality:

```cpp
BoundingBox bbox(atoms, padding);

// Get dimensions and properties
Point3D min = bbox.getMin();
Point3D max = bbox.getMax();
Point3D center = bbox.getCenter();
Point3D dimensions = bbox.getDimensions();
double volume = bbox.getVolume();
double surfaceArea = bbox.getSurfaceArea();

// Get all 8 corners
auto corners = bbox.getCorners();

// Check if point is inside
bool inside = bbox.contains(Point3D(1.0, 2.0, 3.0));
```

### Octree

Efficient 3D spatial partitioning data structure:

```cpp
// Create octree from domain coordinates
Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);

// Create octree from existing bounding box
Octree octree(boundingBox);

// Subdivide octree with various termination conditions
octree.subdivide(maxDepth = 8, minCellSize = 0.001);

// Subdivide with custom occupancy check
octree.subdivide(maxDepth, minCellSize, [](const OctreeNode& node) {
    return node.contains(Point3D(0.5, 0.5, 0.5));
});

// Access root node
const OctreeNode& root = octree.getRoot();
std::cout << "Root center: (" << root.center.x << ", " 
          << root.center.y << ", " << root.center.z << ")" << std::endl;

// Get tree statistics
size_t totalNodes = octree.getNodeCount();
size_t leafNodes = octree.getLeafCount();

// Find containing leaf for a point
const OctreeNode* leaf = octree.findLeaf(Point3D(0.25, 0.25, 0.25));

// Print tree structure
octree.printTree(maxDepth = 3);  // Show up to depth 3
```

#### OctreeNode Properties

Each octree node contains:

```cpp
Point3D min, max;        // Bounding box coordinates
Point3D center;          // Center point
Point3D halfSize;        // Half-size in each dimension
int depth;               // Depth level (root = 0)
bool isLeaf;             // True if no children
std::array<std::unique_ptr<OctreeNode>, 8> children;  // 8 octant children

// Node methods
bool contains(const Point3D& point);  // Check if point is inside
double getVolume();                   // Get node volume
```

### VoxelGrid

Uniform voxelization for mesh generation:

```cpp
// Create voxel grid with atoms
VoxelGrid voxelGrid(atoms, voxelSize, padding);

// Create with existing bounding box
VoxelGrid voxelGrid(boundingBox, atoms, voxelSize);

// Get grid information
double size = voxelGrid.getVoxelSize();
const auto& dims = voxelGrid.getDimensions();  // [nx, ny, nz]
int totalVoxels = voxelGrid.getTotalVoxelCount();
int occupiedVoxels = voxelGrid.getOccupiedVoxelCount();
int emptyVoxels = voxelGrid.getEmptyVoxelCount();

// Access voxels
const auto& occupied = voxelGrid.getOccupiedVoxels();
const auto& empty = voxelGrid.getEmptyVoxels();

// Get specific voxel by grid indices
const Voxel* voxel = voxelGrid.getVoxel(i, j, k);

// Print statistics
voxelGrid.printStatistics();
```

### Hexahedral Mesh Generation

Generate finite element meshes from octrees or voxel grids:

```cpp
#include "biomesh2/OctreeMeshGenerator.hpp"
#include "biomesh2/VoxelMeshGenerator.hpp"

// From octree (adaptive mesh)
Octree octree(boundingBox);
octree.subdivide(3);
HexMesh octreeMesh = OctreeMeshGenerator::generateHexMesh(octree);

// From voxel grid (uniform mesh)
VoxelGrid voxelGrid(atoms, 1.0, 2.0);
HexMesh voxelMesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);

// Access mesh data
std::cout << "Nodes: " << mesh.getNodeCount() << "\n";
std::cout << "Elements: " << mesh.getElementCount() << "\n";

// Access individual nodes and elements
const Point3D& node = mesh.nodes[0];
const std::array<int, 8>& element = mesh.elements[0];
```

#### HexMesh Structure

```cpp
struct HexMesh {
    std::vector<Point3D> nodes;                  // Unique node coordinates
    std::vector<std::array<int, 8>> elements;    // Element connectivity
    
    size_t getNodeCount() const;
    size_t getElementCount() const;
};
```

Each element uses standard hexahedral node ordering (0-3: bottom face, 4-7: top face).

## Mesh Generation Approaches

BioMesh2 provides two complementary approaches for hexahedral mesh generation:

### Octree-Based (Adaptive)
- **Use case**: Adaptive FEM, hierarchical methods
- **Advantages**: Lower element count, locally refined
- **Characteristics**: Non-uniform elements, complex connectivity

### Voxelization-Based (Uniform)
- **Use case**: Uniform analysis, image processing, structured grids
- **Advantages**: Simple implementation, predictable structure, high node sharing
- **Characteristics**: Uniform cubic elements, regular connectivity

See [docs/VoxelizationMeshGeneration.md](docs/VoxelizationMeshGeneration.md) and [docs/OctreeMeshGeneration.md](docs/OctreeMeshGeneration.md) for detailed documentation.

## Testing

The project includes comprehensive unit tests (57 tests) covering:

- Atom constructor functionality
- Atomic specification database operations
- Property assignment correctness
- Bounding box calculations
- PDB parsing accuracy
- Octree spatial operations
- Voxel grid generation and occupancy
- Hexahedral mesh generation

Run tests with:
```bash
cd build
./biomesh2_tests
# or
make test
```

## Example Output

```
BioMesh2 C++ Module - PDB Parser and Bounding Box Calculator
===========================================================

=== Atom Information ===
    ID Element      X (Å)      Y (Å)      Z (Å) Radius (Å)   Mass (Da)
--------------------------------------------------------------------------
     0       N      20.154      16.967      10.000       1.550      14.007
     1       C      19.030      16.200       9.500       1.700      12.011
     2       C      18.500      15.300      10.600       1.700      12.011
     3       O      17.400      14.800      10.500       1.520      15.999
     4       C      17.900      17.100       8.900       1.700      12.011

=== Bounding Box Information ===
Min corner: (14.380, 6.380, 5.700)
Max corner: (25.020, 20.300, 17.200)
Center:     (19.700, 13.340, 11.450)
Dimensions: 10.640 × 13.920 × 11.500 Å
Volume:     1703.251 Ų
Surface:    861.098 Ų
```

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.