#!/bin/bash
if [[ -z "$1" ]]; then
    echo "No action supplied"
elif [[ "$1" == "to" ]]; then
  echo "Upload source code to" $2
  if [[ "$2" == "lanl" ]]; then
    scp -r Makefile src qsub apps lanl:gr-fe.lanl.gov:/usr/projects/cint/cint_sces/ExactDiagonalization/
  else
    rsync -avzh \
    --include-from 'code-include.txt' \
    --exclude-from 'code-exclude.txt' \
    $(pwd) $2:~/GitRepo/
  fi
elif [[ "$1" == "from" ]]; then
  echo "Get data from" $2
  if [[ "$2" == "lanl" ]]; then
    scp -r lanl:gr-fe.lanl.gov:/net/scratch3/chengyanlai/$3 /Users/chengyanlai/data/
  else
    rsync -avzh --update \
    --exclude-from 'data-exclude.txt' \
    $2:~/GitRepo/ExactDiagonalization/data/$3 /Users/chengyanlai/data/
  fi
fi
