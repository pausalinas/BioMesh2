# BioMesh2 - C++ PDB Parser and Molecular Bounding Box Calculator

BioMesh2 is a modern C++ module for parsing PDB (Protein Data Bank) files, enriching atoms with physical properties, computing molecular bounding boxes, and performing efficient 3D spatial operations using octree data structures.

**Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.**

## Working Effectively

### Bootstrap, Build, and Test the Repository:

Run these commands in exact order:

```bash
# Install dependencies (Ubuntu/Debian)
sudo apt update && sudo apt install -y build-essential cmake libgtest-dev

# Create build directory
mkdir build && cd build

# Configure with CMake (takes ~1 second)
cmake ..

# Build the project (takes ~5-7 seconds, NEVER CANCEL)
make -j4

# Run tests (takes <1 second, 35 tests, NEVER CANCEL)
make test

# Run examples to verify functionality
./biomesh2_example
./octree_demo
./biomesh2_example ../data/test_peptide.pdb 2.0
```

**TIMING EXPECTATIONS:**
- **CMake Configuration**: 1-2 seconds. NEVER CANCEL - Set timeout to 60+ seconds.
- **Build Process**: 5-7 seconds fresh build. NEVER CANCEL - Set timeout to 300+ seconds.
- **Test Suite**: <1 second for all 35 tests. NEVER CANCEL - Set timeout to 120+ seconds.
- **Examples**: <1 second each. Set timeout to 60+ seconds.

**CRITICAL BUILD NOTE:**
- If you encounter build errors in `Octree.cpp` related to `OCTREE_CHILD_COUNT`, the code has duplicate loop structures. Fix by ensuring only one loop from `i < 8` exists in the `printTreeRecursive` method.

### Prerequisites

**Required:**
- C++17 compatible compiler (GCC 13.3+, Clang, or MSVC)
- CMake 3.14 or later
- Build tools (`build-essential` on Ubuntu)

**Optional but Recommended:**
- GoogleTest for running tests (`libgtest-dev` on Ubuntu)

### Test Validation

**Always run the complete test suite after making changes:**
- Test suite covers: Atom, AtomicSpec, AtomBuilder, BoundingBox, PDBParser, Octree
- 35 comprehensive tests covering all functionality
- All tests must pass for the build to be considered successful

**Manual validation scenarios:**
1. **PDB Parsing**: `./biomesh2_example ../data/test_peptide.pdb 1.5` - should parse 15 atoms
2. **Built-in Demo**: `./biomesh2_example` - should show step-by-step demonstration with 5 atoms
3. **Octree Functionality**: `./octree_demo` - should create 585 nodes with 512 leaves at depth 3
4. **Bounding Box Calculation**: Verify output includes volume, surface area, and correct dimensions

## Key Projects and Components

### Core Library (`biomesh2_lib`)
- **Location**: `src/` and `include/biomesh2/`
- **Purpose**: Core functionality for PDB parsing and spatial operations
- **Components**:
  - `Atom.hpp/cpp` - Enhanced atom class with physical properties
  - `AtomicSpec.hpp/cpp` - Atomic specifications database (H, C, N, O, P, S, F, Cl, Br, I, Na, Mg, K, Ca, Fe, Zn, Se)
  - `PDBParser.hpp/cpp` - PDB file parsing with element extraction
  - `AtomBuilder.hpp/cpp` - Atom property enrichment
  - `BoundingBox.hpp/cpp` - 3D bounding box calculations
  - `Octree.hpp/cpp` - Spatial partitioning data structure
  - `BioMesh2.hpp/cpp` - Main convenience API

### Example Programs
- **`biomesh2_example`**: Demonstrates PDB parsing, atom enrichment, and bounding box calculation
- **`octree_demo`**: Shows octree subdivision and spatial queries

### Test Suite (`biomesh2_tests`)
- **Location**: `tests/test_biomesh2.cpp`
- **Coverage**: All components with 35 comprehensive tests
- **Runtime**: <1 second for complete test suite

## Sample Data
- **Location**: `data/test_peptide.pdb`
- **Content**: Small peptide with 15 atoms (N, C, O elements)
- **Usage**: Perfect for testing and validation

## Common Development Tasks

### Building After Code Changes
```bash
cd build
make -j4  # Takes 2-7 seconds depending on changes
```

### Testing After Changes
```bash
cd build
make test        # Run via CTest
./biomesh2_tests  # Run directly (shows detailed output)
```

### Adding New Chemical Elements
1. Edit `src/AtomicSpec.cpp` in `initializeDefaultSpecs()` method
2. Add new `AtomicSpec` with element symbol, radius, and mass
3. Rebuild and test: element should be recognized in PDB parsing

### Working with PDB Files
- **Supported Format**: Standard PDB ATOM records
- **Element Detection**: Uses columns 77-78 or fallback to atom name parsing
- **Coordinates**: Extracts X, Y, Z from standard PDB columns
- **Example Usage**: `auto atoms = PDBParser::parsePDBFile("file.pdb");`

## API Quick Reference

### One-Line Usage
```cpp
auto [atoms, boundingBox] = processPDBFile("protein.pdb", 2.0); // 2.0 Å padding
```

### Step-by-Step Usage
```cpp
// Parse PDB
auto basicAtoms = PDBParser::parsePDBFile("protein.pdb");

// Enrich with properties
AtomBuilder builder;
auto enrichedAtoms = builder.buildAtoms(basicAtoms);

// Calculate bounding box
BoundingBox bbox(enrichedAtoms, 1.5); // 1.5 Å padding
```

### Octree Operations
```cpp
// Create from bounding box
Octree octree(boundingBox);

// Subdivide with termination conditions
octree.subdivide(maxDepth = 8, minCellSize = 0.001);

// Query operations
const OctreeNode* leaf = octree.findLeaf(Point3D(x, y, z));
```

## Troubleshooting

### Build Issues
- **Missing GoogleTest**: Build succeeds but tests aren't available. Install `libgtest-dev` and reconfigure.
- **Octree.cpp Errors**: Check for duplicate loop code or undefined constants. Ensure only standard C++ constructs.
- **C++17 Errors**: Verify compiler supports C++17. Use GCC 7+ or Clang 5+.

### Runtime Issues
- **PDB Parsing Fails**: Check file format and permissions. Use sample file to verify functionality.
- **Unknown Element Error**: Element not in atomic database. Add to `AtomicSpec.cpp` or use supported elements.
- **Empty Bounding Box**: Ensure atoms have valid coordinates and elements are supported.

## Expected Output Examples

### biomesh2_example Output
```
=== Atom Information ===
    ID Element      X (Å)      Y (Å)      Z (Å) Radius (Å)   Mass (Da)
--------------------------------------------------------------------------
     0       N      20.154      16.967      10.000       0.560      14.007
     1       C      19.030      16.200       9.500       0.670      12.011
...

=== Bounding Box Information ===
Min corner: (14.920, 6.920, 6.230)
Max corner: (24.480, 19.770, 16.670)
Volume:     1282.512 Ų
Surface:    713.613 Ų
```

### Test Output
```
[==========] Running 35 tests from 6 test suites.
...
[  PASSED  ] 35 tests.
```

**Always validate your changes by running through at least one complete end-to-end scenario after making modifications.**