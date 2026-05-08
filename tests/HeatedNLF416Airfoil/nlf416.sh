#!/bin/bash

# Generate UCD file for NLF416 airfoil from a coarse set of airfoil coordinates

export SPAN=1.0
export SPANCELLS=30
export SPANCOORD=y
export AOA=10.19

# Generate GMSH .geo file:
python3 nlf416.py $SPAN $SPANCELLS $SPANCOORD $AOA

# Generate GMSH .msh file:
gmsh -2 -format msh22 -log gmsh.log nlf416.geo

# Convert GMSH .msh file to UCD file:
python3 gmsh2ucd.py nlf416 

mv nlf416.inp ./IBM/nlf416
