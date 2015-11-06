#! /bin/bash

mpif77 *.f -c *.o
ar rcs libSPARSKIT.a *.o
cp libSPARSKIT.a $BIOTISSUE_DIR/install/lib
