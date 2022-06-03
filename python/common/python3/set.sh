#!/bin/bash
#
# Usage: source set.sh
# Python3 must be in your path.

top=`pwd`
alias cg_spica="python3 $top/src/cg_spica.py"
source src/spica-completion.sh

## Directory that is in users' paths.
#bin_dir=$HOME/bin
## Name of script used to excute cg_spica tools.
#install_name=cg_spica
#
#if [ ! -e $bin_dir ];then
#  mkdir $bin_dir
#fi
#chmod u+x src/cg_spica.py
#ln -s $top/src/cg_spica.py $bin_dir/$install_name
