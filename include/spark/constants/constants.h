#pragma once

// source: https://physics.nist.gov/cuu/Constants/Table/allascii.txt

namespace spark::constants {

// Elementary charge [C]
constexpr double e = 1.602176634e-19;

// Electron mass [kg]
constexpr double m_e = 9.1093837139e-31;

// Pi
constexpr double pi = 3.1415926535897932;

// Vacuum permittivity [F m^-1]
constexpr double eps0 = 8.8541878188e-12;

// Boltzmann constant [J K^-1]
constexpr double kb = 1.380649e-23;

// Atomic mass constant [kg]
constexpr double amc = 1.66053906892e-27;

// Helium atomic mass [kg]
constexpr double m_he = 4.002602 * amc;

// Xenon atomic mass [kg]
constexpr double m_xe = 131.293 * amc;

}  // namespace spark::constants