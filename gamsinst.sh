#!/bin/bash
# Copyright (C) 2011 GAMS Development and others
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Michael Bussieck, Stefan Vigerske

gamspath="${1:-gams}"

if test -z "$gamspath" ; then
  echo "Need to get path of GAMS system as first argument."
  exit 1
fi

# TODO adapt for Windows
#gmscmp=${gamspath}/gmscmpNT.txt
gmscmp="${gamspath}/gmscmpun.txt"

gmscmporig="${gmscmp}.orig"
gmscmpbak="${gmscmp}.bak"

bin=`realpath gamsworhp`

if ! test -r "$gmscmp" ;
then
   echo "File $gmscmp not found or not readable, cannot edit."
   exit 1
fi

if ! test -e "$bin"
then
   echo "Solver binary $bin not found, cannot install."
   exit 1
fi

echo "Adding or updating entry for WORHP into $gmscmp (pointing to $bin)."

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
   # TODO adapt for Windows
   #scriptCmd  = "gmsgennt.cmd";
   #execCmd    = "gmsgennx.exe";
   scriptCmd = "gmsworus.run";
   execCmd   = "gmsworux.out";

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
' $gmscmpbak > $gmscmp

cat > ${gamspath}/gmsworus.run <<EOF
#!/bin/sh
${bin} \$4
[ $? = 0 ] || echo "error: solver rc $? != 0" 1>&2; exit 11
EOF
chmod +x ${gamspath}/gmsworus.run

rm -f ${gamspath}/gmsworux.out
ln -s $bin ${gamspath}/gmsworux.out

#echo "Installing $lib in $gamspath"

#if test -e "${gamspath}/$libname" ;
#then
#   rm -f "${gamspath}/$libname"
#fi
#cp "${libdir}/$libname" "${gamspath}/$libname"
#ln -s "${libdir}/$libname" "${gamspath}/$libname"

# hide libstdc++ and some other GCC libraries if system files are more recent than what GAMS ships
if [ -e "${gamspath}/libstdc++.so.6" ] ; then
  # get highest GLIBCXX dependency from libstdc++.so.6 in GAMS
  gamsver=`strings ${gamspath}/libstdc++.so.6 | grep "^GLIBCXX_[0-9]" | sort -V | tail -1`
  # check which libstdc++ is used by gamsworhp
  userlibstdcxx=`ldd gamsworhp | grep libstdc++ | awk '{print $3}'`
  if [ -n "$userlibstdcxx" ] ; then
    # get highest GLIBCXX dependency from libstdc++.so.6 that is used by WORHP
    userver=`strings ${userlibstdcxx} | grep "^GLIBCXX_[0-9]" | sort -V | tail -1`
  fi
  if [ -n "${gamsver}" ] && [ -n "${userver}" ] ; then
    # if things worked so far, compare the version numbers (by checking what sort -V lists last)
    gamsver=${gamsver/GLIBCXX_/}
    userver=${userver/GLIBCXX_/}
    echo "libstdc++ versions GAMS: $gamsver   ${userlibstdcxx}: $userver"
    if [ `printf "$gamsver\n$userver" | sort -V | tail -1` != "$gamsver" ] ; then
      # hide GCC libraries in GAMS system directory
      echo "GAMS libstdc++ is older than dependency of gamsworhp."
      for f in libgcc_s.so.1 libstdc++.so.6 libgfortran.so.3 libquadmath.so.0 ; do
        if [ -e "${gamspath}/$f" ] ; then
          echo "Moving ${gamspath}/$f to ${gamspath}/${f}.hide"
          mv "${gamspath}/$f" "${gamspath}/${f}.hide"
        fi
      done
    else
      echo "GAMS libstdc++ is at least as recent as dependency of gamsworhp, so no action needed."
    fi

  else
    echo "Could not check GLIBCXX version of libstdc++.so.6 in GAMS directory and WORHP library dependency."
    echo "If running GAMS/WORHP fails due some missing GLIBCXX version, try removing libstdc++.so.6 from the GAMS system directory."
  fi

fi
