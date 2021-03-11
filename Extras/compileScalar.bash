#!/usr/bin/bash

# compile the source
gcc $1.c -lm -o $1.exe

# make it executable
chmod 700 $1.exe
