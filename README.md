# elphanharmonie

This project serves two purposes:

* Be faster than `elphmod.md.Driver`
* Help me learn C89

## Installation

The following command compiles the program using GCC without optimization:

    make

Different compilers and optimization flags can be selected:

    make CC=icx CFLAGS=-O2

LAPACK and BLAS are required.

## Input file

The program `elphy` accepts a single argument with the path to an input file
with the following content, where `[…]` are placeholders for numerical values.
The indices i, j, k, l run over lattice vectors, α, β, γ, δ over orbitals, and
x, y, z over the three Cartesian displacement directions for all atoms.

    [temperature kT]
    [number of electrons per unit cell]
    [number of orbitals per unit cell]
    [supercell size]
    [number of directions per unit cell]
    [number of lattice vectors]
    [R₀₀] [R₀₁] [R₀₂]
    [R₁₀] [R₁₁] [R₁₂] ← lattice vectors in units of primitive vectors
    [R₂₀] [R₂₁] [R₂₂]
    ⋮
    [number of hopping parameters]
    [i₀] [α₀] [β₀] [<0 α₀|H|Rᵢ₀ β₀>]
    [i₁] [α₁] [β₁] [<0 α₁|H|Rᵢ₁ β₁>] ← hopping parameters
    [i₂] [α₂] [β₂] [<0 α₂|H|Rᵢ₂ β₂>]
    ⋮
    [number of interatomic force constants]
    [j₀] [x₀] [y₀] [∂²E/[∂u(0, x₀) ∂u(Rⱼ₀, y₀)]
    [j₁] [x₁] [y₁] [∂²E/[∂u(0, x₁) ∂u(Rⱼ₁, y₁)] ← force constants
    [j₂] [x₂] [y₂] [∂²E/[∂u(0, x₂) ∂u(Rⱼ₂, y₂)]
    ⋮
    [number of electron-phonon matrix elements]
    [k₀] [z₀] [l₀] [γ₀] [δ₀] [<0 γ₀|∂H/∂u(Rₖ₀, z₀)|Rₗ₀ δ₀>]
    [k₁] [z₁] [l₁] [γ₁] [δ₁] [<0 γ₁|∂H/∂u(Rₖ₁, z₁)|Rₗ₁ δ₁>] ← matrix elements
    [k₂] [z₂] [l₂] [γ₂] [δ₂] [<0 γ₂|∂H/∂u(Rₖ₂, z₂)|Rₗ₂ δ₂>]
    ⋮
