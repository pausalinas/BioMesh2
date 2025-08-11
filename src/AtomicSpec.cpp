#include "biomesh2/AtomicSpec.hpp"
#include <stdexcept>

namespace biomesh2 {

AtomicSpecDatabase& AtomicSpecDatabase::getInstance() {
    static AtomicSpecDatabase instance;
    return instance;
}

AtomicSpecDatabase::AtomicSpecDatabase() {
    initializeDefaultSpecs();
}

void AtomicSpecDatabase::initializeDefaultSpecs() {
    // Atomic radii (not van der Waals) and atomic masses for common elements
    // Atomic radii data sources: Slater's rules and empirical crystallographic data
    // Atomic masses from NIST atomic masses (2020)
    
    specs_["H"] = AtomicSpec("H", 0.31, 1.008);    // Hydrogen
    specs_["C"] = AtomicSpec("C", 0.67, 12.011);   // Carbon
    specs_["N"] = AtomicSpec("N", 0.56, 14.007);   // Nitrogen
    specs_["O"] = AtomicSpec("O", 0.48, 15.999);   // Oxygen
    specs_["P"] = AtomicSpec("P", 0.98, 30.974);   // Phosphorus
    specs_["S"] = AtomicSpec("S", 0.88, 32.065);   // Sulfur
    specs_["F"] = AtomicSpec("F", 0.42, 18.998);   // Fluorine
    specs_["Cl"] = AtomicSpec("Cl", 0.79, 35.453); // Chlorine
    specs_["Br"] = AtomicSpec("Br", 0.94, 79.904); // Bromine
    specs_["I"] = AtomicSpec("I", 1.15, 126.904);  // Iodine
    specs_["Na"] = AtomicSpec("Na", 1.54, 22.990); // Sodium
    specs_["Mg"] = AtomicSpec("Mg", 1.30, 24.305); // Magnesium
    specs_["K"] = AtomicSpec("K", 1.96, 39.098);   // Potassium
    specs_["Ca"] = AtomicSpec("Ca", 1.74, 40.078); // Calcium
    specs_["Fe"] = AtomicSpec("Fe", 1.17, 55.845); // Iron
    specs_["Zn"] = AtomicSpec("Zn", 1.25, 65.38);  // Zinc
    specs_["Se"] = AtomicSpec("Se", 1.03, 78.96);  // Selenium
}

const AtomicSpec& AtomicSpecDatabase::getSpec(const std::string& element) const {
    auto it = specs_.find(element);
    if (it == specs_.end()) {
        throw std::runtime_error("Element '" + element + "' not found in atomic specification database");
    }
    return it->second;
}

bool AtomicSpecDatabase::hasElement(const std::string& element) const {
    return specs_.find(element) != specs_.end();
}

void AtomicSpecDatabase::addSpec(const AtomicSpec& spec) {
    specs_[spec.elementSymbol] = spec;
}

} // namespace biomesh2