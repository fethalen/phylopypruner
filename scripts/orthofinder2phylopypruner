#!/bin/bash
#
## orthofinder2phylopypruner.bash: Convert OrthoFinder MSAs into PhyloPyPruner-friendly MSAs.
##
## Author: Felix Thalen <felix.thalen@uni-goettingen.de>
##
## Usage:
##   ./orthofinder2phylopypruner.bash <species_list> <MSA directory>
##
## About:
## This script will format the output of OrthoFinder into a format that is
## suitable for PhyloPyPruner. 'OTU_identifier' will be converted into
## 'OTU@identifier'.
##
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

readonly ARGS=( "$@" )
readonly NO_OF_ARGS="${#}"
readonly FASTA_EXT=".fa"
readonly NEWICK_EXT=".tre"

# Print usage instructions.
print_usage_instructions() {
  [ "$*" ] && echo "$0: $*"
  sed -n '/^##/,/^$/s/^## \{0,1\}//p' "$0"
  exit 2
} 2>/dev/null

# Print usage instructions and exit if the number of arguments is not 2.
parse_args() {
  if [[ ${NO_OF_ARGS} != 2 || "${ARGS[0]}" =~ -?-h(elp)? ]]
  then
    print_usage_instructions
    exit 1
  fi
}

# Takes the path to a directory and a list of species as an input.
rename_alignments() {
  local species_list=$1
  local target_dir=$2

  find "${target_dir}" -type f -name "*${FASTA_EXT}" \
    | while read -r filename
  do
    while read -r otu
    do
      sed -i "s/>${otu}_/>${otu}@/g" "${filename}"
    done < "${species_list}"
  done
}

rename_trees() {
  local species_list=$1
  local target_dir=$2

  find "${target_dir}" -type f -name "*${NEWICK_EXT}" \
    | while read -r filename
  do
    while read -r otu
    do
      sed -i "s/${otu}_/(${otu}@/g" "${filename}"
    done < "${species_list}"
  done
}

main() {
  parse_args

  local species_list=${ARGS[0]}
  local target_dir=${ARGS[1]}

  rename_alignments ${species_list} ${target_dir}
  rename_trees ${species_list} ${target_dir}
}

main
