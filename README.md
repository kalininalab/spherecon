# SphereCon
This tool calculates a measure for relative solvent accessible area for single residues or whole amino acid chains.
## Dependecies
* python 3
* numpy
* biopython
## Instructions
### Usage
spherecon.py -i /Path/To/Input/File [-o /Path/To/Output/File] [-c chain] [-r residues] [--ca] [--bb]

-i:     Path to an input file in PDB (Protein Data Bank) file format.

-o:     Path to the output file produced as tab separated text file.
        Default: *InputFile*_spherecon.tsv

-c:     Chain identifier of the amino acid chain, which should be analysed denoted as the chain identifiers of the ATOM records in the PDB file.
        Default: The first chain found in the input file

-r:     List of residue identifiers for all residues for which the SphereCon value should be computed. If not given, the SphereCon values for all residues are computed.
        The residue identifiers denote as the residues identifiers of the ATOM records in the PDB file.
        Examples: 234,78,368 | 17 | 34,35,36,37

--ca:   C alpha version of SphereCon. Needs only coordinates of the C alpha atoms and their amino acid type.

--bb:   Backbone only version of SphereCon. Needs only coordinates of the C alpha atoms.
### Using as library
spherecon can be simply imported by
import spherecon
(just add the file to your PYTHONPATH)
