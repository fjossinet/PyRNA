# PyRNA

<img src="logo.png" width="200px">

A Python library modeling and using RNA concepts

## Installation and usage
To use PyRNA in your scripts, clone this repository and add its location to your PYTHONPATH env variable.
To use the PyRNA scripts from any location, add the scripts dir to your PATH env variable.
Something like:
<pre>
export PYRNA_HOME=$HOME/PyRNA
export PYTHONPATH=$PYTHONPATH:$PYRNA_HOME
export PATH=$PATH:$PYRNA_HOME/scripts
</pre>

## PyRNA Library
the source code for PyRNA is divided into 3 modules:
* model: the RNA concepts. So far, the available concepts are:
    * Block: contiguous positions between a start and an end position
    * Location: a list of Block objects
    * Molecule: stores a name, an organism name, a sequence and a bunch of modified residues
    * RNA: a specialized Molecule. Its sequence can be made with the following characters : A, U, G, C, - and _. 
    * BasePair
    * SecondaryStructure
    * Atom: stores an atom name, x, y and z coordinates
    * TertiaryStructure
* parsers: functions to load PyRNA objects from RNA files (PDB, FASTA, Vienna, RNAML,...) or to dump them into files
* db: load data from public databases


## PyRNA scripts
scripts using the PyRNA library

## More to come soon
The package will be further extended according to the educational demand