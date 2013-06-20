#!/bin/sh

#  featureExtractionShell.sh
#
#  Created by Drenkow, Nathan G. on 12/3/12.
#  Copyright (c) 2012 JHU/APL. All rights reserved.

# ========= EXAMPLE ==========
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib 
#PATCH_DIR="/home/Testing/"
#PATCH_FILENAME="universalPatches400PerSize.xml"
#IMAGE_DIR="/home/Testing/images/"
#EXECUTABLE="/home/Testing/HMAX_FE"
# ============================

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib #OPENCV LIBRARY FOLDER LOCATION
export LD_LIBRARY_PATH

PATCH_DIR=""
PATCH_FILENAME=""
IMAGE_DIR=""
EXECUTABLE=""

# Order of argument pairs does not matter
$EXECUTABLE "-id" $IMAGE_DIR "-pd" $PATCH_DIR "-pf" $PATCH_FILENAME
