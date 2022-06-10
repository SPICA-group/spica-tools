#!/bin/bash
#
# Usage: source set.sh
# Python3 must be in your path.
SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
if [ ! -e $SCRIPT_DIR/bin ];then
  mkdir $SCRIPT_DIR/bin
  ln -s $SCRIPT_DIR/src/cg_spica.py $SCRIPT_DIR/bin/cg_spica
fi
export PATH=$PATH:$SCRIPT_DIR/bin
source $SCRIPT_DIR/src/spica-completion.sh
