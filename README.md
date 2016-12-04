Starting point for a GAMS/Worhp interface.

Create a symlink 'gams' pointing to a GAMS system directory.

Create a symlink 'worhp' pointing to a Worhp build.
Directories worhp/lib and worhp/include should exist.

Run make and make install.
The latter will install the build of the GAMS/Worhp interface as solver
WORHP in GAMS. The GAMS system directory need to be writable.
