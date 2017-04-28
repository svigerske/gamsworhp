#!/bin/sh
set -e

GAMSPATH="/path/to/gams"
WORHPPATH="libworhp_vs2013_1.9-1_16_10/vs2013-Release"

GAMSAPI="$GAMSPATH/apifiles/C/api"

cl -nologo -c "$GAMSAPI/gmomcc.c" -I"$GAMSAPI"
cl -nologo -c "$GAMSAPI/gevmcc.c" -I"$GAMSAPI"
#TODO SRCDIR is wrong, but should not be used with Worhp 2 anymore
cl -nologo -c "gamsworhp.c" -I"$GAMSAPI" -I"$WORHPPATH/include/worhp" -DSRCDIR=\"\"

link -nologo -out:gamsworhp.exe gmomcc.obj gevmcc.obj gamsworhp.obj "$WORHPPATH"/bin32/libworhp.lib

./gamsinstwin.sh "$GAMSPATH" "$WORHPPATH"
