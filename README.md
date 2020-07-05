# DerivatizeMe

## Installation

This requires the presence of openbabel-2, and is tested with openbabel 2.4.1. For older openbabel 2.0 you will need to comment out the line:
311     themolecule.SetFormula(formula);//
Edit Makefile accordingly if there are special locations for openbabel 2 include and library locations.
type "make"
for use type "./dzme"

## Usage

```
Correct usage:

           ./dzme -m coremolecule -h hydrogensubs -I substituentdir -O outputdir -Of outputformat -L logfile -s listofsubstituents -x sdfoutput -m maximumsubs

Example:   ./dzme -m ./molecule.xyz -h list 7,9,13 -I ./substituents -O ./outputdir -L ./logfile.txt
            -s cl.xyz,no2.xyz,ch3.xyz -x output.sdf -i 

 -m coremolecule:
                    core molecule to be derivatized (default "core.pdb")
 -h hydrogensubs:
    hydrogensubs:   list 7,9,13     {only allows substitution of hydrogens 7,9,13 systematic}
                                    {note that this must be the exact atom numbers in the entire molecule}
                    all             {allows substitution of all hydrogens in a systematic manner}
                    ring            {allows substitution of all hydrogens attached to rings a systematic manner}
                    random  30      {allows substitution of approximately thirty percent of all hydrogens systematically (default)
                    systematic      {allows monosubstitution of all hydrogens in a systematic manner}
 -I substituentdir:
                    directory containing substituents (default "./subs")
 -O outputdir:
                    directory containing individual files created (default "./output")
 -Of outputformat:
                    file type for individual files created (default "none")
                    set as "none" to only create an sdf file output
 -L logfile:
                    log file containing smiles/energies of all structures
                    default "logfile.txt"
 -s substituents:
                    list of structures in the substituent directory that will be used for functionalization
                    default is "cl.pdb,nh2.pdb"
 -x outputsdf:
                    sdf file containing all generated structures (default "output.sdf")
 -i:
                    for all and ring, only make C-H's substitutable (default)
                    ignored for -h list
 -a:
                    for all and ring, make all H's substitutable
                    ignored for -h list
 -maxs maxsubstituents:
                    the maximum number of substituents in the products (default 100)
 -maxt maxtotal:
                    cut off production when maxtotal is reached (default 100000)
 -opt:
                    perform MMFF optimization(default off)
 -phys:
                    calculate physicochemical properties(default off)
 -noforce:
                    allow sites to remain unsubstituted(default force unless maxs set)
 -conf:
                    search for a better conformer(default off)
```
## Examples


