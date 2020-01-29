# SphereCon
This tool calculates an approximate values for relative solvent accessible area for single residues or whole amino acid chains in a PDB file. The method is based on calculating or estimating the number of atoms around a given amino acid residue that are located within a sphere with a cut-out cone. 
## Dependencies
* python 3
* numpy
* biopython
## Instructions
### Usage
spherecon.py -i /Path/To/Input/File [-o /Path/To/Output/File] [-c chain] [-r residues] [--ca] [--bb] [--dm]
-i:     Path to an input file in PDB (Protein Data Bank) file format or distance matrix file in case of --dm.

-o:     Path to the output file produced as tab separated text file.
        Default: *InputFile*_spherecon.tsv

-c:     Chain identifier of the amino acid chain, or a list of Chain identifier separated by a ',', which should be analysed denoted as the chain identifiers of the ATOM records in the PDB file.
        Default: The first chain found in the input file
        Examples: A | A,B | C,B,z,2

-r:     List of residue identifiers for all residues for which the SphereCon value should be computed. If not given, the SphereCon values for all residues are computed.
        The residue identifiers denote as the residues identifiers of the ATOM records in the PDB file.
        Examples: 234,78,368 | 17 | 34,35,36,37

--ca:   C alpha version of SphereCon. Needs only coordinates of the C alpha atoms and their amino acid type.

--bb:   Backbone only version of SphereCon. Needs only coordinates of the C alpha atoms.

--dm:   Distance matrix version of SphereCon. Takes a distance matrix as input instead of a PDB file.
        File for distance matrix has to be in the following format:
                SEQ [amino acid sequence]
                [residue nr 1] [residue nr 2] [distance value in angstrom]
                ...
        Example:
                SEQ TGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSN
                4 9 12.641
                5 8 10.656
                5 9 12.616
                5 10 13.342
                5 37 11.964

### Using as library
spherecon can be simply imported by
import spherecon
(just add the file to your PYTHONPATH)
## Cite
Gress A, Kalinina OV (2019) SphereCon - A method for precise estimation of residue relative solvent accessible area from limited structural information. Submitted.
