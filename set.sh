#!/bin/bash
#
# Usage: source set.sh
# Python3 must be in your path.
SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
export PATH=$PATH:$SCRIPT_DIR/bin
source $SCRIPT_DIR/src/spica-completion.sh
