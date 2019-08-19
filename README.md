# Membrane-Proteins-Tutorial-Simulation-using-NAMD-and-VMD

#############################################################

Main Source of the tutorial:
http://www.ks.uiuc.edu/Training/Tutorials/science/membrane/mem-tutorial.pdf

IMP --> Warning, this tutorial is designed for users with limited computational resources (2 fs timestep, the SHAKE algorithm for H atoms, and a multiple time-stepping algorithm for the computation of long-range electrostatic interactions).

Ideally, one should use a uniform 1 fs time step without SHAKE. This is particularly important for long simulations performed at constant energy.

However, you may want to keep water molecules rigid, as the CHARMM force field was parameterized using rigid water molecules.

#############################################################

(1) Setting up a structural model of a membrane protein starting from a raw PDB file.

(1.1) Viewing and examining the Protein.

VMD Tcl commands:
```
mol new NPSR1_Asn107Ile_raw.pdb
```
Notes:
- Different structural elements are divided into “chains” in the PDB file. look at the COMPND item.
- Use the proper representation.
- Use the Query mouse mode to see information about individual atoms.
- All structures are missing some atoms and residues . (Imp)
- You should always delete selections once you are finished with them. (Imp)
- VMD uses 4×4 matrices to perform rotations and translations, so one must add a fourth row ({0 0 0 1}) to the matrices provided in the PDB file.

(1.2) Generating PSF and PDB Files (Preparing PDB contains the coordinates alone without hydrogen)

VMD Tcl commands:
```
mol new NPSR1_Asn107Ile_raw.pdb
set NPSR1 [atomselect top protein]
$NPSR1 writepdb NPSR1_Asn107Ile.pdb
```
-create PSF
-load pdb.then, use auto PSF generator and use (top all27 prot lipid.rtf) or C:/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf

VMD Tcl commands:
```
mol new NPSR1_Asn107Ile.pdb
```

(1.3) Water box (Use add solvation box or command line)

VMD Tcl commands:
```
file delete NPSR1_Asn107Ile.pdb
file delete NPSR1_Asn107Ile_autopsf.log
file delete NPSR1_Asn107Ile_autopsf_formatted.pdb

mol new NPSR1_Asn107Ile_autopsf.psf
mol addfile NPSR1_Asn107Ile_autopsf.pdb

solvate NPSR1_Asn107Ile_autopsf.psf NPSR1_Asn107Ile_autopsf.pdb -t 5 -o NPRS1_Sol_raw
```

(1.4) Remove unwanted water molecule (use hydrophobicity (by visualizing ResType)).

- Charged residues are colored blue and red, green residues are hydrophilic, and white residues are hydrophobic.
- Change the representation involving water molecules to cover different regions along the z-axis.

- To make selection easier, we will first center everything on origin. Then, water molecule selection.

VMD Tcl commands:
```
mol new NPRS1_Sol_raw.psf
mol addfile NPRS1_Sol_raw.pdb
set all [atomselect top all]
$all moveby [vecinvert [measure center $all]]
display resetview
```
- Then, create the proper graphical representations (check notes). Now we will “mark” the water we don’t want to keep using thr B-factor field.

VMD Tcl commands:
```
set solv [atomselect top "segname WT1"]
$solv set beta 1
```
- unmark the waters we want to keep.

VMD Tcl commands:
```
set seltext "segname WT1 and same residue as \
((y < -20) or (y > 20))"
set sel [atomselect top $seltext]
$sel set beta 0
```
- Make a list of these for the marked molecules.

VMD Tcl commands:
```
set badwater [atomselect top "name OH2 and beta > 0"]
set seglist [$badwater get segid]
set reslist [$badwater get resid]
```
- Now we use psfgen to read the PSF and PDB files, and delete the residues we just marked:

VMD Tcl commands:
```
mol delete all
package require psfgen
resetpsf
topology C:/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf

readpsf NPRS1_Sol_raw.psf
coordpdb NPRS1_Sol_raw.pdb

foreach segid $seglist resid $reslist {
delatom $segid $resid
}

writepdb NPRS1_Sol_bilayer.pdb
writepsf NPRS1_Sol_bilayer.psf
```
- We now have a complete protein with water molecules in just the right places.

#------------------------------------------------------------------------------------------

(2) The steps needed to place the protein in a native-like membrane environment.

(2.1) Building a Membrane Patch (use Membrane Builder (name-->POPC))  --> x=80 y=80 topology=c27 name=popc

(2.2) Alignment of Membrane and Protein

VMD Tcl commands:
```
set popc [atomselect top all]
set NPRS1mol [mol new NPRS1_Sol_bilayer.psf]
mol addfile NPRS1_Sol_bilayer.pdb
set NPRS1 [atomselect $NPRS1mol all]
$popc moveby [vecinvert [measure center $popc weight mass]]
$popc writepdb popc_TEMP.pdb

set vest [atomselect $NPRS1mol "protein"]
$NPRS1  moveby [vecinvert [measure center $vest weight mass]]
display resetview

$NPRS1 move [transaxis x -90]
$NPRS1 writepdb NPRS1_temp.pdb
```
- combine the two temporary files into one set of PSF and PDB files:

VMD Tcl commands:
```
mol delete all
package require psfgen
resetpsf
readpsf popc.psf
coordpdb popc_TEMP.pdb
readpsf NPRS1_Sol_bilayer.psf
coordpdb NPRS1_temp.pdb

writepsf NPRS1_popc_raw.psf
writepdb NPRS1_popc_raw.pdb

file delete popc_TEMP.pdb
file delete NPRS1_temp.pdb
```
- Load new files to visualize them.

(2.3) Combination of Membrane and Protein -->  Make room for PDB in the membrane layer, so that the protein doesn’t overlap any lipid molecules.

- We will again use the beta field of the atoms to mark the “bad lipids.”

VMD Tcl commands:
```
mol delete all
mol new NPRS1_popc_raw.psf
mol addfile NPRS1_popc_raw.pdb
set POPC "resname POPC"
set all [atomselect top all]
$all set beta 0

set seltext1 "$POPC and same residue as \
(name P1 and z>0 and abs(x)<15 and abs(y)<15)"
set seltext2 "$POPC and same residue as \
(name P1 and z<0 and abs(x)<10 and abs(y)<10)"
set seltext3 "$POPC and same residue as (within 0.6 of protein)"
set sel1 [atomselect top $seltext1]
set sel2 [atomselect top $seltext2]
set sel3 [atomselect top $seltext3]
$sel1 set beta 1
$sel2 set beta 1
$sel3 set beta 1
set badlipid [atomselect top "name P1 and beta > 0"]
set seglistlipid [$badlipid get segid]
set reslistlipid [$badlipid get resid]
```

- We must remove the water molecules overlapping with the protein:

VMD Tcl commands:
```
set seltext4 "(water and not segname WCA WCB WCC WCD WF SOLV) \
and same residue as within 3 of \
((same residue as (name P1 and beta>0)) or protein)"
set seltext5 "segname SOLV and same residue as \
within 3 of lipids"
set sel4 [atomselect top $seltext4]
set sel5 [atomselect top $seltext5]
$sel4 set beta 1
$sel5 set beta 1
set badwater [atomselect top "name OH2 and beta > 0"]
set seglistwater [$badwater get segid]
set reslistwater [$badwater get resid]
```
- Now, we delete the atoms as we did before:

VMD Tcl commands:
```
mol delete all
resetpsf
readpsf NPRS1_popc_raw.psf
coordpdb NPRS1_popc_raw.pdb
foreach segid $seglistlipid resid $reslistlipid {
delatom $segid $resid
}
foreach segid $seglistwater resid $reslistwater {
delatom $segid $resid
}
writepsf NPRS1_popc.psf
writepdb NPRS1_popc.pdb

file delete NPRS1_Sol_bilayer.psf
file delete NPRS1_Sol_bilayer.pdb
file delete NPRS1_popc_raw.psf
file delete NPRS1_popc_raw.pdb
file delete popc.psf
file delete popc.pdb
```
(2.4) Solvation and Ionization of the whole system.

- First, let’s use the minmax option of the measure command to see how big our water layer is:

VMD Tcl commands:
```
mol delete all
mol new NPRS1_popc.psf
mol addfile NPRS1_popc.pdb
set water [atomselect top water]
measure minmax $water
````
- use minmax to guide the dim of the box.
- The -b 1.5 option tells solvate to remove atoms within 1.5 °A of the solute

VMD Tcl commands:
```
package require solvate
solvate NPRS1_popc.psf NPRS1_popc.pdb -o NPRS1_popc_water_TEMP -b 1.5 -s WA \
-minmax {{-50 -50 -50} {50 50 50}}

file delete NPRS1_popc.psf
file delete NPRS1_popc.pdb
```
- Remove water inside the lipid bilayer and around the protein atoms.

VMD Tcl commands:
```
mol delete all
mol new NPRS1_popc_water_TEMP.psf
mol addfile NPRS1_popc_water_TEMP.pdb

set all [atomselect top all]
$all set beta 0
set seltext "segid WT1 to WT99 and same residue as abs(z) < 25"
set sel [atomselect top $seltext]
$sel set beta 1
set badwater [atomselect top "name OH2 and beta > 0"]
set seglist [$badwater get segid]
set reslist [$badwater get resid]

mol delete all
package require psfgen
resetpsf
readpsf NPRS1_popc_water_TEMP.psf
coordpdb NPRS1_popc_water_TEMP.pdb
foreach segid $seglist resid $reslist {
delatom $segid $resid
}
writepdb NPRS1_popcw.pdb
writepsf NPRS1_popcw.psf

file delete NPRS1_popc_water_TEMP.psf
file delete NPRS1_popc_water_TEMP.pdb
file delete NPRS1_popc_water_TEMP.log
```

- Ionization

VMD Tcl commands:
```
mol delete all
mol new NPRS1_popcw.psf
mol addfile NPRS1_popcw.pdb
```

- Extensions -> Modeling -> Add Ions (check the notes)
- Output prefix = NPRS1_popcwi
- salt -> NaCl

VMD Tcl commands:
```
mol new NPRS1_popcwi.psf
mol addfile NPRS1_popcwi.pdb
```
#------------------------------------------------------------------------------------------

(3) Minimize and equilibrate the resulting system with NAMD (4 stages -> Continuing jobs).

(3.1) Melting of Lipid Tails (simulation in which everything except lipid tails, is fixed to induce the appropriate disorder of a fluid-like bilayer )

VMD Tcl commands:
```
mol delete all
mol new NPRS1_popcwi.psf
mol addfile NPRS1_popcwi.pdb
set all [atomselect top "all"]
$all set beta 0
set fixed [atomselect top "water or name CLA POT or protein or \
(chain L and name O2 P1 O3 O4 O1 C15 H52 H51 H11 C11 H12 \
N C14 H42 H43 H41 C12 H22 H23 H21 C13 H33 H31 H32)"]
$fixed set beta 1
$all writepdb NPRS1_popcwi.fix
exit
```

- Stage 1 (check RUN_NPRS1_popcwimineq-01.conf)
- use findbox.tcl -> chnage Periodic Boundary Conditions in RUN_NPRS1_popcwimineq-01.conf.
- check Fixed Atoms Constraint -> name.fix
- Minimization -> 1000 (NAMD will run 1000 steps of minimization, then it will reinitiate velocities according to the desired temperature)
- run 250000 ;# 0.5 ns (using a 2 fs timestep)

VMD Tcl commands:
```
mol new NPRS1_popcwi.psf
mol addfile NPRS1_popcwi.pdb
source findbox.tcl
get_cell
exit
```

```
namd2 +p8 RUN_NPRS1_popcwimineq-01.conf > RUN_NPRS1_popcwimineq-01.log.out
```

(3.2) Minimization and Equilibration with Protein Constrained.

- The system has many unnatural atomic positions.
- second run will be a “minimization” run, which guides the system to the nearest local energy minimum in configuration space.
- Minimization will be then followed by equilibration with the protein constrained, so as to permit the environment to relax first.

- Before performing the simulation, we need to generate a file for the harmonic constraints on the protein.
- Choose (fix watexcept {})

VMD Tcl commands:
```
mol new NPRS1_popcwi.psf
mol addfile NPRS1_popcwi.pdb
set all [atomselect top "all"]
$all set beta 0
set prot [atomselect top "protein"]
$prot set beta 1
$all writepdb NPRS1_popcwi.cnst
exit
```

```
namd2 +p8 RUN_NPRS1_popcwimineq-02.conf > RUN_NPRS1_popcwimineq-02.log.out
```

(3.3) Equilibration with Protein Released.

- After minimization and equilibration with the protein constrained, we hopefully have a system in which lipids are well packed around the protein, while water has not entered forbidden regions -> next -> equilibrate the whole system.
- We have eliminated the minimization step in the present simulation.

VMD Tcl commands:
```
namd2 +p8 RUN_NPRS1_popcwimineq-03.conf > RUN_NPRS1_popcwimineq-03.log.out
```

- After the simulation is done, your system should be equilibrated fairly well. You can monitor the stability of the protein through the computation of RMSDs and by looking at the resulting trajectory with VMD.

(3.4) Production Runs

- Now that the protein has been equilibrated, we are ready to perform production runs.
- There will be one main difference with the previous simulations ->  # Constant Pressure Control (useConstantArea command = yes)

VMD Tcl commands:
```
namd2 +p8 RUN_NPRS1_popcwimineq-04.conf > RUN_NPRS1_popcwimineq-04.log.out
```

Simulation of membrane bilayers has shown that the current CHARMM forcefield parameters do no reproduce the experimentally observed area per lipid over long MD trajectories. During the previous simulation steps, we let the are in the xy-plane fluctuate to permit packing of lipids against the protein. However, after good packing has been observed, one should keep the area in the xy-plane constant. Ideally, one should also compute the effective area per lipid in the system by computing the area of the simulations cell and subtracting from it the protein area.
