#include <gtest/gtest.h>
#include "biomesh2/Atom.hpp"
#include "biomesh2/AtomicSpec.hpp"
#include "biomesh2/AtomBuilder.hpp"
#include "biomesh2/BoundingBox.hpp"
#include "biomesh2/PDBParser.hpp"
#include <memory>
#include <vector>

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
    
    // Test some specific values
    const auto& carbon = db.getSpec("C");
    EXPECT_EQ("C", carbon.elementSymbol);
    EXPECT_EQ(1.70, carbon.radius);
    EXPECT_EQ(12.011, carbon.mass);
    
    const auto& hydrogen = db.getSpec("H");
    EXPECT_EQ("H", hydrogen.elementSymbol);
    EXPECT_EQ(1.20, hydrogen.radius);
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
    
    // Check carbon properties
    const auto& enrichedCarbon = enrichedAtoms[0];
    EXPECT_EQ("C", enrichedCarbon->getChemicalElement());
    EXPECT_EQ(1.0, enrichedCarbon->getX());
    EXPECT_EQ(2.0, enrichedCarbon->getY());
    EXPECT_EQ(3.0, enrichedCarbon->getZ());
    EXPECT_EQ(1.70, enrichedCarbon->getAtomicRadius());
    EXPECT_EQ(12.011, enrichedCarbon->getAtomicMass());
    EXPECT_EQ(1, enrichedCarbon->getId());
    
    // Check nitrogen properties
    const auto& enrichedNitrogen = enrichedAtoms[1];
    EXPECT_EQ("N", enrichedNitrogen->getChemicalElement());
    EXPECT_EQ(4.0, enrichedNitrogen->getX());
    EXPECT_EQ(5.0, enrichedNitrogen->getY());
    EXPECT_EQ(6.0, enrichedNitrogen->getZ());
    EXPECT_EQ(1.55, enrichedNitrogen->getAtomicRadius());
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