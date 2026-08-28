# ![elphy](logo.svg)

`elphy`, short for electron-phonon anharmonicity, is a C program to calculate
free energies and forces for tight-binding models on supercells, complemented
with linear electron-phonon coupling and a harmonic potential. It is a faster
rewrite of the Python module `elphmod.md` and thus primarily a driver for the
path-integral molecular-dynamics code i-PI. It is written in ANSI C (C89/C99).

The source-code repository is found at <https://codeberg.org/janberges/elphy>.

## Installation

The following command compiles the program using GCC without optimization:

    make

Different compilers and optimization flags can be selected:

    make CC=icx CFLAGS=-O2

LAPACK and BLAS with the standard LP64 interface are required.

## Usage

The program accepts one, two, or three arguments:

    elphy <data file>
    elphy <data file> <socket>
    elphy <data file> <number> <radius>

With one argument, it alternately reads atomic positions in the XYZ format from
standard input and writes the supercell vectors, atomic positions, free energy,
and forces in ASE's extended XYZ format to standard output.

With two arguments, it exchanges these quantities with i-PI through its socket
interface. `<socket>` is a host name optionally followed by a colon and a port
number or by `/shm`. If a port number is given, an internet socket is used for
communication. Otherwise, a UNIX socket and, optionally, shared memory is used.
The address must match the information in the i-PI input file `input.xml`.

With three arguments, it prints `<number>` sets of atomic positions with random
displacements smaller than `<radius>` and the corresponding energies and forces
in ASE's extended XYZ format. If `<number>` is negative, only the positions are
output using i-PI's standard XYZ format.

The `<data file>` is defined below:

    <temperature kT>
    <number of electrons per unit cell>
    <number of orbitals per unit cell>
    <number of spins per orbital>
    <strain>
    A₀₀ A₀₁ A₀₂
    A₁₀ A₁₁ A₁₂
    A₂₀ A₂₁ A₂₂
    a₀₀ a₀₁ a₀₂
    a₁₀ a₁₁ a₁₂
    a₂₀ a₂₁ a₂₂
    <number of atoms per unit cell>
    X₀ r₀₀ r₀₁ r₀₂ F₀₀ F₀₁ F₀₂
    X₁ r₁₀ r₁₁ r₁₂ F₁₀ F₁₁ F₁₂
    X₂ r₂₀ r₂₁ r₂₂ F₂₀ F₂₁ F₂₂
    ⋮
    <number of lattice vectors>
    R₀₀ R₀₁ R₀₂
    R₁₀ R₁₁ R₁₂
    R₂₀ R₂₁ R₂₂
    ⋮
    <number of hopping parameters>
    i₀ α₀ β₀ <0 α₀|H|Rᵢ₀ β₀>
    i₁ α₁ β₁ <0 α₁|H|Rᵢ₁ β₁>
    i₂ α₂ β₂ <0 α₂|H|Rᵢ₂ β₂>
    ⋮
    <number of interatomic force constants>
    j₀ x₀ y₀ ∂²E/[∂u(0, x₀) ∂u(Rⱼ₀, y₀)]
    j₁ x₁ y₁ ∂²E/[∂u(0, x₁) ∂u(Rⱼ₁, y₁)]
    j₂ x₂ y₂ ∂²E/[∂u(0, x₂) ∂u(Rⱼ₂, y₂)]
    ⋮
    <number of electron-phonon matrix elements>
    k₀ z₀ l₀ γ₀ δ₀ <0 γ₀|∂H/∂u(Rₖ₀, z₀)|Rₗ₀ δ₀>
    k₁ z₁ l₁ γ₁ δ₁ <0 γ₁|∂H/∂u(Rₖ₁, z₁)|Rₗ₁ δ₁>
    k₂ z₂ l₂ γ₂ δ₂ <0 γ₂|∂H/∂u(Rₖ₂, z₂)|Rₗ₂ δ₂>
    ⋮

The indices `i, j, k, l` run over lattice vectors, `α, β, γ, δ` over orbitals,
and `x, y, z` over the three Cartesian displacement directions for all atoms.
All indices are zero-based. All matrix elements `<…>` are single real numbers.
The primitive, position, and force vectors `a, r, F` are given in Cartesian,
the supercell and lattice vectors `A, R` in integer crystal coordinates.

No unit conversions are performed, so any consistent energy and length units can
be used. However, i-PI expects energies and forces in Hartree atomic units. The
atom labels `X` are copied to the `<init file>` and define the default masses.

Note that the interatomic force constants are assumed to be partially screened:
They shall exclude the harmonic term of the electronic potential-energy surface.
Any forces that the model may generate at zero displacements can be compensated
by adding a force correction specified next to the atomic positions.

Lattice vectors and atomic positions are multiplied by 1 + `<strain>` and the
zero-displacement hopping parameters and lattice energy adjusted accordingly.

## Tests and examples

The makefile provides some recipes that exemplify the usage of the program:

- `make test model=TaS2` (or `graphene`) creates `<data file>` and `<init file>`
  for a model and verifies that the computed free energy and forces are correct.
- `make ipi` lets  `elphy` and i-PI perform a structural relaxation together.
- `make show` displays an animation of the relaxation process.
- `make clean` removes compiled files, `make distclean` all generated files.

The Python packages `elphmod`, `ipi`, and `matplotlib` are required.

## Acknowledgment

The method is described in the following paper:

- Arne Schobert, Jan Berges, Erik G. C. P. van Loon, Michael A. Sentef, Sergey
  Brener, Mariana Rossi, and Tim O. Wehling, [SciPost Phys. **16**, 046 (2024)
  ](https://doi.org/10.21468/SciPostPhys.16.2.046)
