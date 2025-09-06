#include <gtest/gtest.h>
#include "biomesh2/Atom.hpp"
#include "biomesh2/AtomicSpec.hpp"
#include "biomesh2/AtomBuilder.hpp"
#include "biomesh2/BoundingBox.hpp"
#include "biomesh2/PDBParser.hpp"
#include "biomesh2/Octree.hpp"
#include "biomesh2/OctreeMeshGenerator.hpp"
#include <memory>
#include <vector>
#include <fstream>
#include <sstream>

using namespace biomesh2;

// Test fixtures
class AtomTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class AtomicSpecTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class AtomBuilderTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class BoundingBoxTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class PDBParserTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class OctreeTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class OctreeMeshGeneratorTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Atom class tests
TEST_F(AtomTest, ConstructorWithElementOnly) {
    Atom atom("C");
    EXPECT_EQ("C", atom.getChemicalElement());
    EXPECT_EQ(0.0, atom.getX());
    EXPECT_EQ(0.0, atom.getY());
    EXPECT_EQ(0.0, atom.getZ());
    EXPECT_EQ(0.0, atom.getAtomicRadius());
    EXPECT_EQ(0.0, atom.getAtomicMass());
    EXPECT_EQ(0, atom.getId());
}

TEST_F(AtomTest, ConstructorWithElementAndRadius) {
    Atom atom("N", 1.55);
    EXPECT_EQ("N", atom.getChemicalElement());
    EXPECT_EQ(1.55, atom.getAtomicRadius());
    EXPECT_EQ(0.0, atom.getAtomicMass());
}

TEST_F(AtomTest, ConstructorWithAllProperties) {
    Atom atom("O", 1.52, 15.999);
    EXPECT_EQ("O", atom.getChemicalElement());
    EXPECT_EQ(1.52, atom.getAtomicRadius());
    EXPECT_EQ(15.999, atom.getAtomicMass());
}

TEST_F(AtomTest, SetCoordinates) {
    Atom atom("H");
    atom.setCoordinates(1.0, 2.0, 3.0);
    EXPECT_EQ(1.0, atom.getX());
    EXPECT_EQ(2.0, atom.getY());
    EXPECT_EQ(3.0, atom.getZ());
}

TEST_F(AtomTest, SetId) {
    Atom atom("C");
    atom.setId(42);
    EXPECT_EQ(42, atom.getId());
}

// AtomicSpec tests
TEST_F(AtomicSpecTest, DatabaseSingleton) {
    auto& db1 = AtomicSpecDatabase::getInstance();
    auto& db2 = AtomicSpecDatabase::getInstance();
    EXPECT_EQ(&db1, &db2);
}

TEST_F(AtomicSpecTest, DefaultSpecsExist) {
    auto& db = AtomicSpecDatabase::getInstance();
    
    // Test common elements
    EXPECT_TRUE(db.hasElement("H"));
    EXPECT_TRUE(db.hasElement("C"));
    EXPECT_TRUE(db.hasElement("N"));
    EXPECT_TRUE(db.hasElement("O"));
    EXPECT_TRUE(db.hasElement("P"));
    EXPECT_TRUE(db.hasElement("S"));
    
    // Test some specific values - using atomic radii (not van der Waals)
    const auto& carbon = db.getSpec("C");
    EXPECT_EQ("C", carbon.elementSymbol);
    EXPECT_EQ(0.67, carbon.radius);  // Atomic radius
    EXPECT_EQ(12.011, carbon.mass);
    
    const auto& hydrogen = db.getSpec("H");
    EXPECT_EQ("H", hydrogen.elementSymbol);
    EXPECT_EQ(0.31, hydrogen.radius);  // Atomic radius
    EXPECT_EQ(1.008, hydrogen.mass);
}

TEST_F(AtomicSpecTest, UnknownElementThrows) {
    auto& db = AtomicSpecDatabase::getInstance();
    EXPECT_FALSE(db.hasElement("Xx"));
    EXPECT_THROW(db.getSpec("Xx"), std::runtime_error);
}

TEST_F(AtomicSpecTest, AddCustomSpec) {
    auto& db = AtomicSpecDatabase::getInstance();
    AtomicSpec custom("X", 2.0, 100.0);
    db.addSpec(custom);
    
    EXPECT_TRUE(db.hasElement("X"));
    const auto& retrieved = db.getSpec("X");
    EXPECT_EQ("X", retrieved.elementSymbol);
    EXPECT_EQ(2.0, retrieved.radius);
    EXPECT_EQ(100.0, retrieved.mass);
}

// AtomBuilder tests - focusing on property assignment
TEST_F(AtomBuilderTest, CorrectPropertyAssignment) {
    // Create basic atoms
    std::vector<std::unique_ptr<Atom>> basicAtoms;
    auto carbon = std::make_unique<Atom>("C");
    carbon->setCoordinates(1.0, 2.0, 3.0);
    carbon->setId(1);
    
    auto nitrogen = std::make_unique<Atom>("N");
    nitrogen->setCoordinates(4.0, 5.0, 6.0);
    nitrogen->setId(2);
    
    basicAtoms.push_back(std::move(carbon));
    basicAtoms.push_back(std::move(nitrogen));
    
    // Build enriched atoms
    AtomBuilder builder;
    auto enrichedAtoms = builder.buildAtoms(basicAtoms);
    
    ASSERT_EQ(2, enrichedAtoms.size());
    
    // Check carbon properties - using atomic radii
    const auto& enrichedCarbon = enrichedAtoms[0];
    EXPECT_EQ("C", enrichedCarbon->getChemicalElement());
    EXPECT_EQ(1.0, enrichedCarbon->getX());
    EXPECT_EQ(2.0, enrichedCarbon->getY());
    EXPECT_EQ(3.0, enrichedCarbon->getZ());
    EXPECT_EQ(0.67, enrichedCarbon->getAtomicRadius());  // Atomic radius
    EXPECT_EQ(12.011, enrichedCarbon->getAtomicMass());
    EXPECT_EQ(1, enrichedCarbon->getId());
    
    // Check nitrogen properties - using atomic radii
    const auto& enrichedNitrogen = enrichedAtoms[1];
    EXPECT_EQ("N", enrichedNitrogen->getChemicalElement());
    EXPECT_EQ(4.0, enrichedNitrogen->getX());
    EXPECT_EQ(5.0, enrichedNitrogen->getY());
    EXPECT_EQ(6.0, enrichedNitrogen->getZ());
    EXPECT_EQ(0.56, enrichedNitrogen->getAtomicRadius());  // Atomic radius
    EXPECT_EQ(14.007, enrichedNitrogen->getAtomicMass());
    EXPECT_EQ(2, enrichedNitrogen->getId());
}

TEST_F(AtomBuilderTest, UnsupportedElementThrows) {
    std::vector<std::unique_ptr<Atom>> basicAtoms;
    auto unknownAtom = std::make_unique<Atom>("Xx");
    basicAtoms.push_back(std::move(unknownAtom));
    
    AtomBuilder builder;
    EXPECT_THROW(builder.buildAtoms(basicAtoms), std::runtime_error);
}

TEST_F(AtomBuilderTest, UnsupportedElementDetection) {
    std::vector<std::unique_ptr<Atom>> basicAtoms;
    basicAtoms.push_back(std::make_unique<Atom>("C"));
    basicAtoms.push_back(std::make_unique<Atom>("Xx"));
    basicAtoms.push_back(std::make_unique<Atom>("Yy"));
    
    AtomBuilder builder;
    EXPECT_FALSE(builder.areAllElementsSupported(basicAtoms));
    
    auto unsupported = builder.getUnsupportedElements(basicAtoms);
    EXPECT_EQ(2, unsupported.size());
    EXPECT_TRUE(std::find(unsupported.begin(), unsupported.end(), "Xx") != unsupported.end());
    EXPECT_TRUE(std::find(unsupported.begin(), unsupported.end(), "Yy") != unsupported.end());
}

// BoundingBox tests - focusing on correct calculation
TEST_F(BoundingBoxTest, CorrectBoundingBoxCalculation) {
    // Create test atoms with known coordinates and radii
    std::vector<std::unique_ptr<Atom>> atoms;
    
    // Atom at origin with radius 1.0
    auto atom1 = std::make_unique<Atom>("C", 1.0, 12.0);
    atom1->setCoordinates(0.0, 0.0, 0.0);
    atoms.push_back(std::move(atom1));
    
    // Atom at (10, 0, 0) with radius 2.0
    auto atom2 = std::make_unique<Atom>("N", 2.0, 14.0);
    atom2->setCoordinates(10.0, 0.0, 0.0);
    atoms.push_back(std::move(atom2));
    
    // No padding
    BoundingBox bbox(atoms, 0.0);
    
    // Expected bounds:
    // Atom1 at (0,0,0) with radius 1.0: bounds (-1,-1,-1) to (1,1,1)
    // Atom2 at (10,0,0) with radius 2.0: bounds (8,-2,-2) to (12,2,2)
    // Combined: Min: (-1, -2, -2), Max: (12, 2, 2)
    Point3D expectedMin(-1.0, -2.0, -2.0);
    Point3D expectedMax(12.0, 2.0, 2.0);
    
    EXPECT_NEAR(expectedMin.x, bbox.getMin().x, 1e-6);
    EXPECT_NEAR(expectedMin.y, bbox.getMin().y, 1e-6);
    EXPECT_NEAR(expectedMin.z, bbox.getMin().z, 1e-6);
    EXPECT_NEAR(expectedMax.x, bbox.getMax().x, 1e-6);
    EXPECT_NEAR(expectedMax.y, bbox.getMax().y, 1e-6);
    EXPECT_NEAR(expectedMax.z, bbox.getMax().z, 1e-6);
}

TEST_F(BoundingBoxTest, BoundingBoxWithPadding) {
    std::vector<std::unique_ptr<Atom>> atoms;
    auto atom = std::make_unique<Atom>("C", 1.0, 12.0);
    atom->setCoordinates(0.0, 0.0, 0.0);
    atoms.push_back(std::move(atom));
    
    double padding = 2.0;
    BoundingBox bbox(atoms, padding);
    
    // Expected bounds with padding:
    // Min: (-1-2, -1-2, -1-2) = (-3, -3, -3)
    // Max: (1+2, 1+2, 1+2) = (3, 3, 3)
    EXPECT_NEAR(-3.0, bbox.getMin().x, 1e-6);
    EXPECT_NEAR(-3.0, bbox.getMin().y, 1e-6);
    EXPECT_NEAR(-3.0, bbox.getMin().z, 1e-6);
    EXPECT_NEAR(3.0, bbox.getMax().x, 1e-6);
    EXPECT_NEAR(3.0, bbox.getMax().y, 1e-6);
    EXPECT_NEAR(3.0, bbox.getMax().z, 1e-6);
}

TEST_F(BoundingBoxTest, BoundingBoxProperties) {
    std::vector<std::unique_ptr<Atom>> atoms;
    auto atom = std::make_unique<Atom>("C", 0.0, 12.0); // No radius for simple calculation
    atom->setCoordinates(1.0, 2.0, 3.0);
    atoms.push_back(std::move(atom));
    
    BoundingBox bbox(atoms, 1.0); // 1.0 padding
    
    // Expected bounds: (0, 1, 2) to (2, 3, 4)
    Point3D expectedCenter(1.0, 2.0, 3.0);
    Point3D expectedDimensions(2.0, 2.0, 2.0);
    double expectedVolume = 8.0;
    double expectedSurfaceArea = 24.0;
    
    Point3D center = bbox.getCenter();
    Point3D dimensions = bbox.getDimensions();
    
    EXPECT_NEAR(expectedCenter.x, center.x, 1e-6);
    EXPECT_NEAR(expectedCenter.y, center.y, 1e-6);
    EXPECT_NEAR(expectedCenter.z, center.z, 1e-6);
    EXPECT_NEAR(expectedDimensions.x, dimensions.x, 1e-6);
    EXPECT_NEAR(expectedDimensions.y, dimensions.y, 1e-6);
    EXPECT_NEAR(expectedDimensions.z, dimensions.z, 1e-6);
    EXPECT_NEAR(expectedVolume, bbox.getVolume(), 1e-6);
    EXPECT_NEAR(expectedSurfaceArea, bbox.getSurfaceArea(), 1e-6);
}

TEST_F(BoundingBoxTest, EmptyAtomVectorThrows) {
    std::vector<std::unique_ptr<Atom>> emptyAtoms;
    EXPECT_THROW(BoundingBox bbox(emptyAtoms), std::runtime_error);
}

// PDB Parser tests
TEST_F(PDBParserTest, ParseSimplePDBContent) {
    std::string pdbContent = 
        "HEADER    TEST                                    01-JAN-70   TEST            \n"
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N  \n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00           C  \n"
        "END                                                                             \n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(2, atoms.size());
    
    // Check first atom
    EXPECT_EQ("N", atoms[0]->getChemicalElement());
    EXPECT_NEAR(20.154, atoms[0]->getX(), 1e-3);
    EXPECT_NEAR(16.967, atoms[0]->getY(), 1e-3);
    EXPECT_NEAR(10.000, atoms[0]->getZ(), 1e-3);
    EXPECT_EQ(0, atoms[0]->getId());
    
    // Check second atom
    EXPECT_EQ("C", atoms[1]->getChemicalElement());
    EXPECT_NEAR(19.030, atoms[1]->getX(), 1e-3);
    EXPECT_NEAR(16.200, atoms[1]->getY(), 1e-3);
    EXPECT_NEAR(9.500, atoms[1]->getZ(), 1e-3);
    EXPECT_EQ(1, atoms[1]->getId());
}

TEST_F(PDBParserTest, EmptyContentThrows) {
    std::string emptyContent = "HEADER    TEST\nEND\n";
    EXPECT_THROW(PDBParser::parsePDBContent(emptyContent), std::runtime_error);
}

TEST_F(PDBParserTest, NonexistentFileThrows) {
    EXPECT_THROW(PDBParser::parsePDBFile("nonexistent_file.pdb"), std::runtime_error);
}

// Enhanced PDB Parser tests for improved element extraction
TEST_F(PDBParserTest, ElementExtractionFromColumns77_78) {
    // Test element extraction when columns 77-78 are present and valid
    std::string pdbContent = 
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N  \n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00           C  \n"
        "ATOM      3  FE  HEM A   2      18.000  15.000   8.000  1.00 30.00          Fe  \n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(3, atoms.size());
    EXPECT_EQ("N", atoms[0]->getChemicalElement());
    EXPECT_EQ("C", atoms[1]->getChemicalElement());
    EXPECT_EQ("Fe", atoms[2]->getChemicalElement());
}

TEST_F(PDBParserTest, ElementExtractionMissingColumns77_78) {
    // Test element extraction when columns 77-78 are missing (shorter lines)
    std::string pdbContent = 
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00\n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00\n"
        "ATOM      3  O   ALA A   1      17.500  14.500   7.000  1.00 20.00\n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(3, atoms.size());
    EXPECT_EQ("N", atoms[0]->getChemicalElement());
    EXPECT_EQ("C", atoms[1]->getChemicalElement());  // CA -> C (alpha carbon in amino acid)
    EXPECT_EQ("O", atoms[2]->getChemicalElement());
}

TEST_F(PDBParserTest, AmbiguousCAResolution) {
    // Test CA disambiguation: CA in amino acid vs Ca (calcium)
    std::string pdbContentAminoAcid = 
        "ATOM      1  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00\n";
    
    std::string pdbContentCalcium = 
        "ATOM      1  CA  CAL A   1      19.030  16.200   9.500  1.00 20.00\n";
    
    auto atomsAA = PDBParser::parsePDBContent(pdbContentAminoAcid);
    auto atomsCa = PDBParser::parsePDBContent(pdbContentCalcium);
    
    ASSERT_EQ(1, atomsAA.size());
    ASSERT_EQ(1, atomsCa.size());
    EXPECT_EQ("C", atomsAA[0]->getChemicalElement());   // Alpha carbon in alanine
    EXPECT_EQ("Ca", atomsCa[0]->getChemicalElement());  // Calcium in non-amino acid
}

TEST_F(PDBParserTest, TwoLetterElementExtraction) {
    // Test extraction of two-letter elements
    std::string pdbContent = 
        "ATOM      1  MG  HEM A   1      20.000  16.000  10.000  1.00 20.00\n"
        "ATOM      2  ZN  ZNC A   2      19.000  15.000   9.000  1.00 20.00\n"
        "ATOM      3  CL  CLA A   3      18.000  14.000   8.000  1.00 20.00\n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(3, atoms.size());
    EXPECT_EQ("Mg", atoms[0]->getChemicalElement());
    EXPECT_EQ("Zn", atoms[1]->getChemicalElement());
    EXPECT_EQ("Cl", atoms[2]->getChemicalElement());
}

TEST_F(PDBParserTest, InvalidElementFallback) {
    // Test fallback when columns 77-78 contain invalid elements
    std::string pdbContent = 
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00          XX  \n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00          YY  \n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(2, atoms.size());
    EXPECT_EQ("N", atoms[0]->getChemicalElement());   // Fallback to atom name
    EXPECT_EQ("C", atoms[1]->getChemicalElement());   // Fallback to atom name (CA -> C in amino acid)
}

TEST_F(PDBParserTest, ElementValidationAgainstDatabase) {
    // Test that only elements in the atomic database are accepted
    // Use 'Q' which is not a real element and shouldn't be in the database
    std::string pdbContent = 
        "ATOM      1  QQ  UNK A   1      20.000  16.000  10.000  1.00 20.00\n";
    
    // This should throw because no valid atoms can be parsed (QQ -> Q, and Q is not in our database)
    EXPECT_THROW(PDBParser::parsePDBContent(pdbContent), std::runtime_error);
}

TEST_F(PDBParserTest, MixedValidInvalidElements) {
    // Test with mixed valid and invalid elements
    // Use 'Q' which is not a real element and shouldn't be in the database
    std::string pdbContent = 
        "ATOM      1  C   ALA A   1      20.000  16.000  10.000  1.00 20.00\n"
        "ATOM      2  QQ  UNK A   2      21.000  17.000  11.000  1.00 20.00\n";
    
    // Should parse only the valid atom (carbon)
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    ASSERT_EQ(1, atoms.size());
    EXPECT_EQ("C", atoms[0]->getChemicalElement());
}

// Octree tests
TEST_F(OctreeTest, RootCellCreation) {
    // Test creating root cell from domain coordinates
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    
    const OctreeNode& root = octree.getRoot();
    EXPECT_EQ(0, root.depth);
    EXPECT_TRUE(root.isLeaf);
    EXPECT_NEAR(0.0, root.min.x, 1e-6);
    EXPECT_NEAR(0.0, root.min.y, 1e-6);
    EXPECT_NEAR(0.0, root.min.z, 1e-6);
    EXPECT_NEAR(1.0, root.max.x, 1e-6);
    EXPECT_NEAR(1.0, root.max.y, 1e-6);
    EXPECT_NEAR(1.0, root.max.z, 1e-6);
    EXPECT_NEAR(0.5, root.center.x, 1e-6);
    EXPECT_NEAR(0.5, root.center.y, 1e-6);
    EXPECT_NEAR(0.5, root.center.z, 1e-6);
    EXPECT_NEAR(0.5, root.halfSize.x, 1e-6);
    EXPECT_NEAR(0.5, root.halfSize.y, 1e-6);
    EXPECT_NEAR(0.5, root.halfSize.z, 1e-6);
}

TEST_F(OctreeTest, RootCellFromBoundingBox) {
    // Create atoms for bounding box
    std::vector<std::unique_ptr<Atom>> atoms;
    auto atom1 = std::make_unique<Atom>("C", 1.0, 12.0);
    atom1->setCoordinates(0.0, 0.0, 0.0);
    atoms.push_back(std::move(atom1));
    
    auto atom2 = std::make_unique<Atom>("N", 1.0, 14.0);
    atom2->setCoordinates(2.0, 2.0, 2.0);
    atoms.push_back(std::move(atom2));
    
    BoundingBox bbox(atoms, 0.0);
    Octree octree(bbox);
    
    const OctreeNode& root = octree.getRoot();
    EXPECT_EQ(0, root.depth);
    EXPECT_TRUE(root.isLeaf);
    EXPECT_NEAR(bbox.getMin().x, root.min.x, 1e-6);
    EXPECT_NEAR(bbox.getMin().y, root.min.y, 1e-6);
    EXPECT_NEAR(bbox.getMin().z, root.min.z, 1e-6);
    EXPECT_NEAR(bbox.getMax().x, root.max.x, 1e-6);
    EXPECT_NEAR(bbox.getMax().y, root.max.y, 1e-6);
    EXPECT_NEAR(bbox.getMax().z, root.max.z, 1e-6);
}

TEST_F(OctreeTest, BasicSubdivision) {
    // Test subdivision to depth 1
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    octree.subdivide(1); // Subdivide to depth 1
    
    const OctreeNode& root = octree.getRoot();
    EXPECT_FALSE(root.isLeaf);
    EXPECT_EQ(9, octree.getNodeCount()); // 1 root + 8 children
    EXPECT_EQ(8, octree.getLeafCount()); // 8 leaf nodes
    
    // Check that all 8 children exist and are at depth 1
    for (size_t i = 0; i < 8; ++i) {
        ASSERT_NE(nullptr, root.children[i]);
        EXPECT_EQ(1, root.children[i]->depth);
        EXPECT_TRUE(root.children[i]->isLeaf);
        EXPECT_NEAR(0.25, root.children[i]->halfSize.x, 1e-6);
        EXPECT_NEAR(0.25, root.children[i]->halfSize.y, 1e-6);
        EXPECT_NEAR(0.25, root.children[i]->halfSize.z, 1e-6);
    }
}

TEST_F(OctreeTest, DeepSubdivision) {
    // Test subdivision to depth 3 as specified in requirements
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    octree.subdivide(3); // Subdivide to depth 3
    
    const OctreeNode& root = octree.getRoot();
    EXPECT_FALSE(root.isLeaf);
    
    // At depth 3: 1 + 8 + 64 + 512 = 585 total nodes
    // Leaves: 8^3 = 512
    EXPECT_EQ(585, octree.getNodeCount());
    EXPECT_EQ(512, octree.getLeafCount());
}

TEST_F(OctreeTest, ContainsPoint) {
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    
    const OctreeNode& root = octree.getRoot();
    
    // Test points inside domain
    EXPECT_TRUE(root.contains(Point3D(0.5, 0.5, 0.5)));
    EXPECT_TRUE(root.contains(Point3D(0.0, 0.0, 0.0))); // boundary
    EXPECT_TRUE(root.contains(Point3D(1.0, 1.0, 1.0))); // boundary
    
    // Test points outside domain
    EXPECT_FALSE(root.contains(Point3D(-0.1, 0.5, 0.5)));
    EXPECT_FALSE(root.contains(Point3D(1.1, 0.5, 0.5)));
    EXPECT_FALSE(root.contains(Point3D(0.5, 0.5, 1.1)));
}

TEST_F(OctreeTest, FindLeaf) {
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    octree.subdivide(2); // Subdivide to depth 2
    
    // Test finding leaf for center point
    Point3D centerPoint(0.5, 0.5, 0.5);
    const OctreeNode* leaf = octree.findLeaf(centerPoint);
    ASSERT_NE(nullptr, leaf);
    EXPECT_TRUE(leaf->isLeaf);
    EXPECT_TRUE(leaf->contains(centerPoint));
    
    // Test finding leaf for corner point
    Point3D cornerPoint(0.0, 0.0, 0.0);
    leaf = octree.findLeaf(cornerPoint);
    ASSERT_NE(nullptr, leaf);
    EXPECT_TRUE(leaf->isLeaf);
    EXPECT_TRUE(leaf->contains(cornerPoint));
    
    // Test point outside domain
    Point3D outsidePoint(-0.1, 0.5, 0.5);
    leaf = octree.findLeaf(outsidePoint);
    EXPECT_EQ(nullptr, leaf);
}

TEST_F(OctreeTest, NodeVolume) {
    Octree octree(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    
    const OctreeNode& root = octree.getRoot();
    EXPECT_NEAR(8.0, root.getVolume(), 1e-6); // 2×2×2 = 8
    
    // Subdivide and check child volumes
    octree.subdivide(1);
    for (size_t i = 0; i < 8; ++i) {
        ASSERT_NE(nullptr, root.children[i]);
        EXPECT_NEAR(1.0, root.children[i]->getVolume(), 1e-6); // Each child: 1×1×1 = 1
    }
}

TEST_F(OctreeTest, MinCellSizeTermination) {
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    
    // Set min cell size to 0.3, which should prevent subdivision beyond depth 1
    // (depth 1 cells have size 0.5, depth 2 cells would have size 0.25 < 0.3)
    octree.subdivide(10, 0.3); // High max depth, but limited by min cell size
    
    EXPECT_EQ(9, octree.getNodeCount()); // 1 root + 8 children (no further subdivision)
    EXPECT_EQ(8, octree.getLeafCount());
}

TEST_F(OctreeTest, OccupancyCheck) {
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    
    // Define occupancy check that only allows subdivision for nodes containing point (0.5, 0.5, 0.5)
    auto occupancyCheck = [](const OctreeNode& node) -> bool {
        return node.contains(Point3D(0.5, 0.5, 0.5));
    };
    
    octree.subdivide(3, 0.001, occupancyCheck);
    
    // Should have fewer nodes than full subdivision since occupancy check limits subdivision
    EXPECT_LT(octree.getNodeCount(), 585);
    EXPECT_GT(octree.getNodeCount(), 1); // But more than just root
}

// OctreeMeshGenerator tests
TEST_F(OctreeMeshGeneratorTest, BasicMeshGeneration) {
    // Create a simple octree with one subdivision
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    octree.subdivide(1); // Creates 8 leaf cells
    
    // Generate mesh
    HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
    
    // Should have 8 elements (one per leaf)
    EXPECT_EQ(8, mesh.getElementCount());
    
    // Should have 27 unique nodes (3x3x3 grid)
    EXPECT_EQ(27, mesh.getNodeCount());
    
    // Verify each element has 8 nodes
    for (const auto& element : mesh.elements) {
        EXPECT_EQ(8, element.size());
        
        // All node indices should be valid
        for (int nodeIndex : element) {
            EXPECT_GE(nodeIndex, 0);
            EXPECT_LT(nodeIndex, static_cast<int>(mesh.getNodeCount()));
        }
    }
}

TEST_F(OctreeMeshGeneratorTest, GiDExportSingleElement) {
    // Create octree with single element
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
    
    // Should have 1 element and 8 nodes
    EXPECT_EQ(1, mesh.getElementCount());
    EXPECT_EQ(8, mesh.getNodeCount());
    
    // Export to temporary file
    std::string filename = "/tmp/test_single_element.msh";
    
    // Test the export function doesn't throw
    EXPECT_NO_THROW(OctreeMeshGenerator::exportToGiD(mesh, filename));
    
    // Verify file was created and has correct content
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());
    
    std::string line;
    
    // Check header
    std::getline(file, line);
    EXPECT_EQ("MESH dimension 3 Elemtype Hexahedra Nnode 8", line);
    
    // Check coordinates section start
    std::getline(file, line);
    EXPECT_EQ("coordinates", line);
    
    // Read 8 coordinate lines
    std::vector<std::string> coordLines;
    for (int i = 0; i < 8; ++i) {
        std::getline(file, line);
        coordLines.push_back(line);
    }
    
    // Check coordinates section end
    std::getline(file, line);
    EXPECT_EQ("end coordinates", line);
    
    // Check elements section start
    std::getline(file, line);
    EXPECT_EQ("elements", line);
    
    // Read element line
    std::getline(file, line);
    // Should start with "1 " (1-based indexing) and have 8 more numbers
    EXPECT_TRUE(line.find("1 ") == 0);
    
    // Count numbers in element line
    std::istringstream iss(line);
    std::string token;
    int count = 0;
    while (iss >> token) {
        count++;
    }
    EXPECT_EQ(9, count); // Element ID + 8 node indices
    
    // Check elements section end
    std::getline(file, line);
    EXPECT_EQ("end elements", line);
    
    file.close();
}

TEST_F(OctreeMeshGeneratorTest, GiDExportMultipleElements) {
    // Create octree with multiple elements
    Octree octree(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    octree.subdivide(1); // Creates 8 leaf cells
    HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
    
    // Export to temporary file
    std::string filename = "/tmp/test_multiple_elements.msh";
    
    // Test the export function doesn't throw
    EXPECT_NO_THROW(OctreeMeshGenerator::exportToGiD(mesh, filename));
    
    // Verify file was created
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());
    
    // Read and parse the file line by line
    std::string line;
    size_t coordCount = 0;
    size_t elemCount = 0;
    bool inCoordinates = false;
    bool inElements = false;
    
    while (std::getline(file, line)) {
        if (line == "coordinates") {
            inCoordinates = true;
            inElements = false;
        } else if (line == "end coordinates") {
            inCoordinates = false;
        } else if (line == "elements") {
            inElements = true;
            inCoordinates = false;
        } else if (line == "end elements") {
            inElements = false;
        } else if (inCoordinates && !line.empty()) {
            // Count non-empty lines in coordinates section
            if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
                coordCount++;
            }
        } else if (inElements && !line.empty()) {
            // Count non-empty lines in elements section
            if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
                elemCount++;
            }
        }
    }
    file.close();
    
    // Verify counts match expected values
    EXPECT_EQ(mesh.getNodeCount(), coordCount);
    EXPECT_EQ(mesh.getElementCount(), elemCount);
}

TEST_F(OctreeMeshGeneratorTest, GiDExportInvalidFilename) {
    // Create a simple mesh
    Octree octree(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    HexMesh mesh = OctreeMeshGenerator::generateHexMesh(octree);
    
    // Try to export to an invalid path
    std::string invalidFilename = "/invalid_directory/test.msh";
    
    // Should throw an exception
    EXPECT_THROW(OctreeMeshGenerator::exportToGiD(mesh, invalidFilename), std::runtime_error);
}

TEST_F(OctreeMeshGeneratorTest, GiDExportEmptyMesh) {
    // Create empty mesh
    HexMesh emptyMesh;
    
    // Export should work even with empty mesh
    std::string filename = "/tmp/test_empty_mesh.msh";
    EXPECT_NO_THROW(OctreeMeshGenerator::exportToGiD(emptyMesh, filename));
    
    // Verify file was created and has correct structure
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());
    
    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());
    file.close();
    
    // Should still have the proper structure
    EXPECT_TRUE(content.find("MESH dimension 3 Elemtype Hexahedra Nnode 8") != std::string::npos);
    EXPECT_TRUE(content.find("coordinates") != std::string::npos);
    EXPECT_TRUE(content.find("end coordinates") != std::string::npos);
    EXPECT_TRUE(content.find("elements") != std::string::npos);
    EXPECT_TRUE(content.find("end elements") != std::string::npos);
}