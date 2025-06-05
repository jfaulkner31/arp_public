#!/bin/bash

###
# Quick script for making documentation so that I do not forget how.
# NOTE: Need sphinx package to compile.
#       pip install sphinx
###
echo
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo NOW BUILDING DOCUMENTATION
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo
cd documentation
make html

cd ../docs
