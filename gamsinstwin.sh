#!/bin/bash
# Copyright (C) 2011 GAMS Development and others
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Michael Bussieck, Stefan Vigerske

gamspath="$1"

if test -z "$gamspath" ; then
  echo "Need to get path of GAMS system as first argument."
  exit 1
fi

gmscmp=${gamspath}/gmscmpNT.txt
gmscmporig="${gmscmp}.orig"
gmscmpbak="${gmscmp}.bak"

if ! test -r "$gmscmp" ;
then
   echo "File $gmscmp not found or not readable, cannot edit."
   exit 1
fi

if ! test -e "gamsworhp.exe"
then
   echo "Solver binary gamsworhp.exe not found, cannot install."
   exit 1
fi

echo "Adding or updating entry for WORHP into $gmscmp."

# keep backup of original gmscmpun.txt file
if ! test -r "$gmscmporig"
then
   cp "$gmscmp" "$gmscmporig"
fi

# keep backup of current gmscmpun.txt file
cp -f "$gmscmp" "$gmscmpbak"

awk '
BEGIN {
   fileType      = 111;
   dictType      = 5;
   licCodes      = "0001020304";
   defaultOkFlag = 1;
   hiddenFlag    = 0;
   scriptCmd = "gmswornt.cmd";
   execCmd   = "gmswornx.exe";

   written["WORHP"] = 0;
   libid["WORHP"] = "wor";
   dicttype["WORHP"] = 5;
   modeltypes["WORHP"] = "LP RMIP QCP RMIQCP NLP DNLP RMINLP CNS";

   startBlock = 0;
}

function writeConfig(solverID) {
   print solverID, fileType, dicttype[solverID], licCodes, defaultOkFlag, hiddenFlag, "1", modeltypes[solverID];
   print scriptCmd;
   print execCmd;
   written[solverID] = 1;
}

(/^\*/ || /^ *$/) { print $0 }

/^DEFAULTS/ {
   for( solverID in written )
      if( !written[solverID] )
      {
         writeConfig(solverID)
         print "";
      }
   print;
   next;
}

!(/^*/ || /^ *$/) {
   if( startBlock < 0 )
   {
      startBlock++;
      next;
   }
   if( $1 in written && !written[$1] )
   {
      writeConfig($1)
      startBlock = -($7+1);
      next;
   }
   print;
}
' "$gmscmpbak" > "$gmscmp"

cat > "${gamspath}/gmswornt.cmd" <<EOF
@echo off
gmswornx.exe "%~4"
if not %errorlevel% == 0 echo err: solver rc %errorlevel% 1>&2
EOF
chmod +x "${gamspath}/gmswornt.cmd"

cp gamsworhp.exe "${gamspath}/gmswornx.exe"
cp $2/bin32/libworhp.dll "${gamspath}"