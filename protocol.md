# [CHARMM][]

[CHARMM]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/index.html

## [Theoretical Background - Molecular Dynamics simulations][theory]

[theory]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node1.html

### [Time evolution][]

[Time evolution]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node8.html

#### [The Verlet algorithm][]

**Questions and answers.**

>   Explain briefly the components of a typical energy function.

There are five components:

*   Bond stretching
    *   Represents covalent bonds between atoms.
    *   Energy is based on the deviation from the ideal distance between atoms.
    *   Usually sufficient to use harmonic approximation.
*   Angle bending
    *   Energy resulting from the angle between any two atoms bound to the same atom.
*   Torsion angles
    *   Rotation of atom groups around a given bond
*   Van-der-Walls terms
    *   Forces not based on covalent or ionic bonds.
    *   Can be approximated by a Lennard-Jones 12-6 potential.
*   Electrostatic term
    *   Non-bonded interactions.
    *   Number of terms is $O(n^2)$ for $n$ atoms.
    *   Cutoff distance can be used to bound computational complexity.

>   Explain the basic ideas underlying the Steepest Descent and Newton Raphson algorithms.

*   Steepest Descent
    *   Start with some coordinates on the energy landscape.
    *   Successively move a discrete amount along the local downhill gradient
        (steepest descent).
    *   Finds a local minimum.
*   Newton
    *   TODO.

>   Explain how the energy function can be used to obtain the dynamics of a protein
>   (Molecular Dynamic Simulation).

*   The force exerted on an any atom can be obtained by taking the derivative of the
    energy function.
*   Newton's second law gives us accelerations.
*   Numerical integration can be used to get approximate trajectories.

>   A typical MD Simulation protocol follows four steps: Initialization, Heating,
>   Equilibration and Production.  What is done in these steps and why is that necessary?

*   Initialization

    TODO

*   Heating

    TODO

*   Equilibration

    TODO

*   Production

    TODO

[The Verlet algorithm]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node9.html

### [Running CHARMM][]

>   Explain briefly which informations are stored in the parameter and topology file.

*   [Parameter file][par_all27_prot_lipid.prm]
    *   From `/sw/sci/app/charmm/c32b1/toppar/par_all27_prot_lipid.prm`.
    *   TODO.
*   [Topology file][top_all27_prot_lipid.rtf]
    *   From `/sw/sci/app/charmm/c32b1/toppar/top_all27_prot_lipid.rtf`.
    *   Contains definitions of biomolecules.

>   Explain briefly what atom types are and why they are needed?

*   Different sets of parameters have to be used to simulate, e.g., a hydrogen atom
    bonded to a carbon atom and one bonded to oxygen.
*   This is called the chemical environment and represented by *atom types*.
*   Each atom can have different types to model its chemical environment.

>   How many atom types are there in an alanine residue?

*   The definition of the alanine residue in the [topology file][top_all27_prot_lipid.rtf]
    is

        ...
        GROUP
        ATOM CB   CT3    -0.27
        ATOM HB1  HA      0.09
        ATOM HB2  HA      0.09
        ATOM HB3  HA      0.09
        GROUP
        ...
    Thus, there are two atom types: `CT3` and `HA`.

>   What is the partial charge of the hydrogen (HS/HG1) atom and the sulfur atom (S/SG)?

*   "First, an atom type is assigned to each atom of the residue. The numerical value
    in the same line defines the partial charge of that atom."  ([1][]).

[1]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node11.html

*   In [top_all27_prot_lipid.rtf][] we find:

        RESI CYS          0.00
        ...
        ATOM SG   S      -0.23
        ATOM HG1  HS      0.16
        ...

    Thus, the partial charges are `0.16` and `-0.23`, respectively.

>   What are the parameters for the sulfur-hydrogen bond (S-HS)?

TODO.

>   What are the parameters for the sulfur-sulfur bond (SM-SM)?

TODO.

>   Which bond is harder to stretch?

TODO.

[Running CHARMM]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node13.html

## [Initialization][]

[Initialization]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node14.html

### [Structure preparation]

**Questions and answers.**

>   Which method was used to determine the structure `1BPI.pdb` experimentally?

[X-ray diffraction][1BPI].

>   Are there other methods to determine experimentally the structure of a protein?

TODO.

>   How many cysteines?  What are their residue numbers?

6 cysteines with residue numbers 5, 14, 30, 38, 51, 55

>   What is the resolution of this structure of BPTI?

[1.09 Å][1BPI].


>   Are there disulfide-bridges in BPTI?

Yes:

    ...
    SSBOND   1 CYS A    5    CYS A   55                          1555   1555  2.01
    SSBOND   2 CYS A   14    CYS A   38                          1555   1555  2.03
    SSBOND   3 CYS A   30    CYS A   51                          1555   1555  2.02
    ...

>   Which temperature was used for the measurements?

[125 K][1BPI].

>   Why are the hydrogen atoms missing?

TODO.  Less visual clutter?

[1BPI]: http://www.rcsb.org/pdb/explore/explore.do?structureId=1BPI

[Structure preparation]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node15.html

[top_all27_prot_lipid.rtf]: top_all27_prot_lipid.rtf
[par_all27_prot_lipid.prm]: par_all27_prot_lipid.prm

<!-- vim: set tw=90 sts=-1 sw=4 et spell: -->
