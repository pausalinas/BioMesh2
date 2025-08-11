#pragma once

#include <string>
#include <unordered_map>

namespace biomesh2 {

/**
 * @brief Atomic specification holding element properties
 */
struct AtomicSpec {
    std::string elementSymbol;
    double radius;  // Atomic radius in Angstroms (not van der Waals radius)
    double mass;    // Atomic mass in Daltons

    AtomicSpec() : elementSymbol(""), radius(0.0), mass(0.0) {}
    AtomicSpec(const std::string& symbol, double r, double m) 
        : elementSymbol(symbol), radius(r), mass(m) {}
};

/**
 * @brief Atomic specifications database
 * 
 * Provides default atomic properties for common elements.
 * Uses atomic radii (not van der Waals radii) for proper atomic representation.
 */
class AtomicSpecDatabase {
public:
    /**
     * @brief Get singleton instance of the database
     */
    static AtomicSpecDatabase& getInstance();

    /**
     * @brief Get atomic specification for an element
     * @param element Chemical element symbol
     * @return AtomicSpec for the element
     * @throws std::runtime_error if element not found
     */
    const AtomicSpec& getSpec(const std::string& element) const;

    /**
     * @brief Check if element exists in database
     * @param element Chemical element symbol
     * @return true if element exists
     */
    bool hasElement(const std::string& element) const;

    /**
     * @brief Add or update atomic specification
     * @param spec AtomicSpec to add/update
     */
    void addSpec(const AtomicSpec& spec);

private:
    AtomicSpecDatabase();
    void initializeDefaultSpecs();

    std::unordered_map<std::string, AtomicSpec> specs_;
};

} // namespace biomesh2