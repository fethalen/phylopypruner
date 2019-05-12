#!/bin/bash
#
## orthofinder2phylopypruner.bash: Convert OrthoFinder MSAs into PhyloPyPruner-friendly MSAs.
##
## Author: Felix Thalen <felix.thalen@uni-goettingen.de>
##
## About:
##   This script will format the output of OrthoFinder into a format that is
##   suitable for PhyloPyPruner. OrthoFinder should be run using the flags '-os -M
##   msa'. Output MSAs from OrthoFinder will be stored in the
##   'Single_Copy_Orthologue_Sequences' directory. Use this directory as your MSA
##   directory. Also provide a list of species names, separated by line, within a
##   file. See usage instructions below.
##
## Usage:
##   ./orthofinder2phylopypruner.bash <species_list> <MSA directory>
##
## What will happen to my files?
##   '>OTU_identifier' will be converted into '>OTU@identifier'
##
## Example <species list> file:
##   Drosophila
##   Caenorhabditis
##   Gallus
##   ...
##
## Example <MSA directory>:
##   OG0000117.fa
##   OG0000118.fa
##   OG0000119.fa
##   ...

set -euo pipefail

# Print instructions.
usage() {
  [ "$*" ] && echo "$0: $*"
  sed -n '/^##/,/^$/s/^## \{0,1\}//p' "$0"
  exit 2
} 2>/dev/null

if [[ $# != 2 ]]; then
  usage "$@"
  exit 1
fi

readonly SPECIES_LIST=$1
readonly TARGET_DIR=$2

find "${TARGET_DIR}" -type f -name '*.fa' | while read -r filename; do
  while read -r otu; do
    gsed -i "s/${otu}_/${otu}@/g" "${filename}"
  done < "${SPECIES_LIST}"
done
