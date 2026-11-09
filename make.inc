DESTDIR = ~/ALADYN/bin

# ANFOMOD, FOCEX
FF = gfortran

# SCOP8, THERMACOND
FC = mpif90
F90 = mpif90

FLAGS = -O2 -fcheck=all -ffree-line-length-0

# the settings configure found on this computer (not in git); they replace those above
-include $(dir $(lastword $(MAKEFILE_LIST)))make.local
