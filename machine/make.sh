#!/bin/sh

# < =========================   DEVELOPER'S NOTE   ========================= > #
#                                                                              #
#  make file   : make.sh                                                       #
#  output file : log_make.txt                                                  #
#                                                                              #
#  This file uses the "make" command to compile corresponding codes            #
#  according to the settings in the "Makefile".                                #
#                                                                              #
#  > make.sh                                                                   #
#                                                                              #
# > ======================================================================== < #


make clean
make -j > log_make.txt
