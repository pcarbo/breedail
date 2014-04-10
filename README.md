#Breeder selection in advanced intercross lines

This repository contains R source code to select breeders in an
advanced intercross line so that inbreeding in the offspring is
minimized.

* [kinship.R](code/kinship.R) defines functions to calculate and
update pairwise kinship coefficients.

* [find_mates.R](code/find_mates.R) pairs breeders in a given
generation of the advanced intercross line in such a way that they
will produce the smallest average inbreeding coefficients in the
offspring.

See [example.R](code/example.R) for an example of how to use these
functions to select mouse breeders in generation F10 of an advanced
intercross line. The pedigree data for this mouse population are found
in [data](data). This script should generate matings similar to those
in [F10pairs.txt](results/F10pairs.txt).


