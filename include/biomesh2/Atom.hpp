#pragma once

#include <string>
#include <memory>

namespace biomesh2 {

/**
 * @brief Enhanced Atom structure with physical properties
 * 
 * Represents an atom with chemical element, coordinates, and physical properties.
 */
class Atom {
public:
    /**
     * @brief Constructor with chemical element only
     * @param element Chemical element symbol (e.g., "C", "N", "O")
     */
    explicit Atom(const std::string& element);

    /**
     * @brief Constructor with chemical element and atomic radius
     * @param element Chemical element symbol
     * @param radius Atomic radius in Angstroms
     */
    Atom(const std::string& element, double radius);

    /**
     * @brief Constructor with chemical element, atomic radius, and atomic mass
     * @param element Chemical element symbol
     * @param radius Atomic radius in Angstroms
     * @param mass Atomic mass in Daltons
     */
    Atom(const std::string& element, double radius, double mass);

    // Getters
    const std::string& getChemicalElement() const { return chemicalElement_; }
    double getX() const { return x_; }
    double getY() const { return y_; }
    double getZ() const { return z_; }
    double getAtomicRadius() const { return atomicRadius_; }
    double getAtomicMass() const { return atomicMass_; }
    size_t getId() const { return id_; }

    // Setters for coordinates
    void setCoordinates(double x, double y, double z);
    void setId(size_t id) { id_ = id; }

private:
    std::string chemicalElement_;
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    double atomicRadius_ = 0.0;
    double atomicMass_ = 0.0;
    size_t id_ = 0;
};

} // namespace biomesh2