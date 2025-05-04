# PyRNA

<img src="logo.png" width="200px">

A Python library modeling and using RNA concepts

## Installation and usage
To use PyRNA in your scripts, clone this repository and add its location to your PYTHONPATH env variable.
To use the PyRNA scripts from any location, add the scripts dir to your PATH  env variable.
Something like:
<pre>
export PYRNA_HOME=$HOME/PyRNA
export PYTHONPATH=$PYTHONPATH:$PYRNA_HOME
export PATH=$PATH:$PYRNA_HOME/scripts
</pre>

## PyRNA Library
the source code for PyRNA is divided into 3 modules:
* model: the RNA concepts
* parsers: functions to load PyRNA objects from RNA files (PDB, FASTA, Vienna, RNAML,...) or to dump them into files
* db: load data from public databases

## PyRNA scripts
scripts using the PyRNA library

