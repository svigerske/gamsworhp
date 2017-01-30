Starting point for a GAMS/Worhp interface.

Create a symlink 'gams' pointing to a GAMS system directory.

Create a symlink 'worhp' pointing to a Worhp build.
Directories worhp/lib and worhp/include should exist.

Run make and make install.
The latter will install the build of the GAMS/Worhp interface as solver
WORHP in GAMS. The GAMS system directory need to be writable.


### Testing

The test/ directory contains a very basic functionality to run GAMS
solvers on a set of instances and produce GAMS trace files, which contain
for each run a line with information on solution time, solve status, etc.

To use it:

1. Change to the test/ directory.
2. Create a symlink called "instances" pointing to a directory with GAMS
   instances, e.g., the "data/gms" directory from MINLPLib 2
   (http://www.gamsworld.org/minlp/minlplib2/minlplib2.zip).
3. Call ```make```. By default, this will run GAMS/WORHP on each instance
   specified in minlplib.test with a timelimit of 30 seconds.

The ```make``` call creates for each solver run on an instance a trace
file in subdirectory trc/, a log file in subdirectory log/, and a GAMS
listing file in subdirectory lst/. The invidual trace files are then
combined into a single one. Remove a trace file to repeat a run.

The following options can be passed to ```make```:
- TESTSET: specify the testset, a file $(TESTSET).test is expected.
- SOLVER: the GAMS solver to be run
- OPTFILE: the solver option file to be used, or 0 for none
- TIME: timelimit
