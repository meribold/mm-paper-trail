% [Molecular modelling][] [practicum][]
  Protocol
% Lukas Waymann

[Molecular modelling]: http://www.bisb.uni-bayreuth.de/Lecture/
[practicum]: http://www.bisb.uni-bayreuth.de/Lecture/practical/practical.html

# CHARMM^[<http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/index.html>]

## Theoretical Background - Molecular Dynamics simulations^[<http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node1.html>]

### Time evolution^[<http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node8.html>]

#### The Verlet algorithm^[<http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node9.html>]

**Questions and answers.**

>   Explain briefly the components of a typical energy function.

There are five components:

*   [Bond stretching](http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node3.html)
    *   Represents covalent bonds between atoms.
    *   Energy is based on the deviation from the ideal distance between atoms.
    *   Usually sufficient to use harmonic approximation.
*   [Angle bending](http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node4.html)
    *   Energy resulting from the angle between any two atoms bound to the same atom.
*   [Torsion angles](http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node5.html)
    *   Rotation of atom groups around a given bond.
*   [Van-der-Walls terms](http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node6.html)
    *   Forces not based on covalent or ionic bonds.
    *   Can be approximated by a Lennard-Jones 12-6 potential.
*   [Electrostatic term](http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node7.html)
    *   Non-bonded interactions.
    *   Number of terms is $O(n^2)$ for $n$ atoms.
    *   Cutoff distance can be used to bound computational complexity.

>   Explain the basic ideas underlying the Steepest Descent and Newton Raphson algorithms.

*   Steepest Descent
    *   Start with some coordinates on the energy landscape.
    *   Successively move a discrete amount along the local downhill gradient (the
        steepest descent).
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

After getting a file describing the structure of a molecule (for example from the
[Protein Data Bank][PDB]), these steps may be performed:

[PDB]: http://www.rcsb.org/pdb

*   Initialization
    *   Add missing hydrogen atoms.
    *   Split the structure into connected segments (partition the [.pdb] file into
        several ones).
    *   Adjust atom names to match those found in the [topology
        file][top_all27_prot_lipid.rtf]
        *   Edit one or the other file.
        *   <http://www.bisb.uni-bayreuth.de/Lecture/Slides/lecture-md.pdf#page=19>
    *   Add water?
    *   Molecule's initial coordinates cause high energies and forces.
    *   

[.pdb]: https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)

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

The definition of the alanine residue in the [topology file][top_all27_prot_lipid.rtf]
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

In the [topology file][top_all27_prot_lipid.rtf] we find:

    RESI CYS          0.00
    ...
    ATOM SG   S      -0.23
    ATOM HG1  HS      0.16
    ...

Knowing that "[t]he numerical value [...] defines the partial charge"[^node11], we see
that the sought after values are `0.16` and `-0.23`, respectively.

[^node11]: <http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node11.html>

>   What are the parameters for the sulfur-hydrogen bond (S-HS)?

    $ grep BONDS -A 7 par_all27_prot_lipid.prm; grep '^S.*HS' par_all27_prot_lipid.prm
    BONDS
    !
    !V(bond) = Kb(b - b0)**2
    !
    !Kb: kcal/mole/A**2
    !b0: A
    !
    !atom type Kb          b0
    S    HS    275.000     1.3250 ! ALLOW   SUL ION

Thus, the force constant $K_b$ for [equation 1][node3] is
$275\ \frac{kcal}{mole\cdot Å^2}$ and $r_{eq}$ for equation 1, which corresponds to `b0`,
is $1.325\ Å$.

[node3]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node3.html
<!--
The relevant lines in the topology file are

    BOND CB HB2  SG HG1
which tells us that a sulfur-hydrogen bond exists, and

    IC CA   CB   SG   HG1   1.5584 113.8700  176.9600  97.1500  1.3341
which tells us the parameters:

*   `1.5584` is the bond length between the carbon atoms.
*   `113.8700` is the angle of the carbon-carbon bond relative to the carbon-sulfur bond.
*   `176.9600` is the torsion angle between the planes defined by the first 3 atoms and
    atoms 2 to 4.

None of those were about the sulfur-hydrogen bond, but the last two are:

*   `97.1500` is the angle between the last 3 atoms, i.e., between the outer carbon-sulfur
    bond and the sulfur-hydrogen bond.
*   `1.3341` is the bond length between the last two atoms: sulfur and hydrogen.
-->

>   What are the parameters for the sulfur-sulfur bond (SM-SM)?

TODO.

>   Which bond is harder to stretch?

TODO.

[Running CHARMM]: http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node13.html

## Initialization^[<http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node14.html>]

### Structure preparation^[<http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node15.html>]

**Questions and answers.**

>   Which method was used to determine the structure `1BPI.pdb` experimentally?

[X-ray diffraction][1BPI].

>   Are there other methods to determine experimentally the structure of a protein?

Other methods include [NMR spectroscopy and electron
microscopy](http://www.bisb.uni-bayreuth.de/Lecture/Slides/lecture-intro.pdf#page=23).

>   How many cysteines?  What are their residue numbers?

6 cysteines with residue numbers 5, 14, 30, 38, 51, 55

>   What is the resolution of this structure of BPTI?

[1.09 Å][1BPI].


>   Are there disulfide-bridges in BPTI?

Yes:

    $ grep '^SSBOND' 1bpi.pdb
    SSBOND   1 CYS A    5    CYS A   55                          1555   1555  2.01
    SSBOND   2 CYS A   14    CYS A   38                          1555   1555  2.03
    SSBOND   3 CYS A   30    CYS A   51                          1555   1555  2.02

>   Which temperature was used for the measurements?

[125 K][1BPI].

>   Why are the hydrogen atoms missing?

<!-- TODO -->
Less visual clutter?  They can be inferred from the incomplete structure?

### Minimization^[<http://www.bisb.uni-bayreuth.de/Lecture/practical/CharmmCourse/Skript/node16.html>]

**Questions and answers.**

>   Why are there high-energy interactions in the crystal structure? 

The molecule may not be in the ideal conformation with respect to minimizing the internal
energy.

>   Why do they have to be removed? 

*   Minima of the energy landscape correspond to stable conformations.
*   Should reduce vibrations in the dynamics simulation.
*   Should get us closer to naturally occuring conformations.
*   <http://www.bisb.uni-bayreuth.de/Lecture/Slides/lecture-mini-nma.pdf#page=16>

>   How much would you expect the atomic coordinates to change during the minimization? 

*   Depends on how close the conformation of the downloaded molecule structure is to an
    energy minimum.
*   Having removed the external water may cause considerable changes.

<!-- http://pandoc.org/MANUAL.html#ending-a-list -->

    $ cp /home/ullmann/Lecture16/PraktMolmod/CharmmCourse/example.inp .
    $ # Copy new parameter and topology files.
    $ # Modify example.inp.
    $ charmm < example.inp > charmm.out
    $ grep 'MINI>' charmm.out > mini.dat
    $ cat mini.dat
    MINI>        0  22814.35962      0.00000    503.70666      0.02000
    MINI>       10    283.03855  22531.32107     94.21325      0.02150
    MINI>       20   -610.28235    893.32089      6.14155      0.00401
    MINI>       30   -685.00044     74.71810      6.92729      0.00431
    MINI>       40   -742.50682     57.50638      1.95547      0.00193
    MINI>       50   -765.12952     22.62270      0.95290      0.00087
    MINI>        0   -765.12952     22.62270      0.95290      0.00000
    MINI>       10   -854.10369     88.97418      1.13469      0.04783
    MINI>       20   -895.41981     41.31612      1.06149      0.04633
    MINI>       30   -911.70797     16.28816      0.68017      0.02692
    MINI>       40   -919.90921      8.20124      0.64534      0.02607
    MINI>       50   -927.86086      7.95165      1.20517      0.03741
    $ head -6 mini.dat > mini-prot.dat
    $ tail -6 mini.dat > mini-watr.dat
    $ xmgrace mini-prot.dat
    $ xmgrace mini-watr.dat

#### Plots after minimization with 50 SD and 50 ABNR steps

##### Protein

![Protein](mini-prot.png "Protein")\ 

##### Water

![Water](mini-watr.png "Water")\ 

[1BPI]: http://www.rcsb.org/pdb/explore/explore.do?structureId=1BPI

[top_all27_prot_lipid.rtf]: top_all27_prot_lipid.rtf
[par_all27_prot_lipid.prm]: par_all27_prot_lipid.prm

<!-- vim: set tw=90 sts=-1 sw=4 et spell: -->