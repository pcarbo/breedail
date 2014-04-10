#Breeder selection in an advanced intercross line

###Overview

This repository contains R functions to select breeders from an
advanced intercross line so that they minimize inbreeding in the
offspring.

* [kinship.R](code/kinship.R) defines functions to calculate and
update kinship coefficients.

* [find_mates.R](code/find_mates.R) pairs breeders in a given
generation of the advanced intercross line that will produce the
smallest average inbreeding coefficients in the offspring.

See [example.R](code/example.R) for an example of how to use these
functions to select mouse breeders in generation F10 of an advanced
intercross line. The pedigree data for this mouse population are found
in the [data](data). This script should generate matings similar to
[F10pairs.txt](results/F10pairs.txt).


