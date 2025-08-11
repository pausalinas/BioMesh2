#include "biomesh2/PDBParser.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace biomesh2 {

std::vector<std::unique_ptr<Atom>> PDBParser::parsePDBFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open PDB file: " + filename);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    return parsePDBContent(buffer.str());
}

std::vector<std::unique_ptr<Atom>> PDBParser::parsePDBContent(const std::string& content) {
    std::vector<std::unique_ptr<Atom>> atoms;
    std::istringstream stream(content);
    std::string line;
    size_t atomId = 0;

    while (std::getline(stream, line)) {
        // Check if line starts with "ATOM"
        if (line.length() >= 4 && line.substr(0, 4) == "ATOM") {
            auto atom = parseAtomLine(line, atomId++);
            if (atom) {
                atoms.push_back(std::move(atom));
            }
        }
    }

    if (atoms.empty()) {
        throw std::runtime_error("No valid ATOM records found in PDB content");
    }

    return atoms;
}

std::unique_ptr<Atom> PDBParser::parseAtomLine(const std::string& line, size_t atomId) {
    // PDB format specification for ATOM records:
    // Columns  Data type      Field        Definition
    // 1-6      Record name    "ATOM  "
    // 7-11     Integer        serial       Atom serial number
    // 13-16    Atom           name         Atom name
    // 17       Character      altLoc       Alternate location indicator
    // 18-20    Residue name   resName      Residue name
    // 22       Character      chainID      Chain identifier
    // 23-26    Integer        resSeq       Residue sequence number
    // 31-38    Real(8.3)      x            Orthogonal coordinates for X in Angstroms
    // 39-46    Real(8.3)      y            Orthogonal coordinates for Y in Angstroms
    // 47-54    Real(8.3)      z            Orthogonal coordinates for Z in Angstroms
    // 77-78    LString(2)     element      Element symbol

    if (line.length() < 54) {
        return nullptr; // Line too short to contain coordinates
    }

    try {
        // Extract coordinates
        double x = parseCoordinate(line, 30, 8);  // positions 31-38 (0-indexed: 30-37)
        double y = parseCoordinate(line, 38, 8);  // positions 39-46 (0-indexed: 38-45)
        double z = parseCoordinate(line, 46, 8);  // positions 47-54 (0-indexed: 46-53)

        // Extract element symbol
        std::string element;
        if (line.length() >= 78) {
            // Try to get element from columns 77-78
            element = line.substr(76, 2);
            // Remove whitespace
            element.erase(std::remove_if(element.begin(), element.end(), ::isspace), element.end());
        }
        
        if (element.empty() && line.length() >= 16) {
            // Fallback: extract from atom name (columns 13-16)
            std::string atomName = line.substr(12, 4);
            element = extractElement(atomName);
        }

        if (element.empty()) {
            return nullptr; // Could not determine element
        }

        // Create atom with just the element (no properties yet)
        auto atom = std::make_unique<Atom>(element);
        atom->setCoordinates(x, y, z);
        atom->setId(atomId);

        return atom;
    } catch (const std::exception&) {
        return nullptr; // Parsing failed
    }
}

std::string PDBParser::extractElement(const std::string& atomName) {
    // Remove whitespace
    std::string name = atomName;
    name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end());
    
    if (name.empty()) {
        return "";
    }

    // For most cases, the first character is the element
    std::string element;
    element += std::toupper(name[0]);
    
    // Check for two-letter elements
    if (name.length() > 1) {
        char second = std::tolower(name[1]);
        // Common two-letter elements in proteins
        if ((element == "C" && second == 'l') ||  // Cl
            (element == "B" && second == 'r') ||  // Br
            (element == "M" && second == 'g') ||  // Mg
            (element == "C" && second == 'a') ||  // Ca
            (element == "F" && second == 'e') ||  // Fe
            (element == "Z" && second == 'n') ||  // Zn
            (element == "S" && second == 'e')) {  // Se
            element += second;
        }
    }

    return element;
}

double PDBParser::parseCoordinate(const std::string& line, size_t start, size_t length) {
    if (start + length > line.length()) {
        throw std::runtime_error("Line too short for coordinate extraction");
    }

    std::string coordStr = line.substr(start, length);
    
    // Remove leading/trailing whitespace
    coordStr.erase(0, coordStr.find_first_not_of(" \t"));
    coordStr.erase(coordStr.find_last_not_of(" \t") + 1);
    
    if (coordStr.empty()) {
        throw std::runtime_error("Empty coordinate field");
    }

    return std::stod(coordStr);
}

} // namespace biomesh2