#!/bin/bash
. /opt/OpenFOAM/OpenFOAM-v2012/etc/bashrc WM_NCOMPROCS=2; export WM_COMPILE_OPTION=Debug
/usr/bin/gdb "$@"