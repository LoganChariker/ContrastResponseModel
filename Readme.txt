0 - Introduction


This package consists of Java source code and Matlab scripts, which together can be used to simulate the macaque V1 model.  All Java code is to be compiled in Java version 7, the latest version of Java compatible with Matlab (the latest version being 2018b at the time of writing).

The primary file containing the model definition is Network.java.  The Java class 'Network' defined there contains fields for model parameters and state, methods for building the model, and methods for evolving model state in time.  A few other helper Java classes are also included.

To build and run the model, we instantiate an object of class Network within the Matlab workspace, and from there fill in model parameters, call model building methods, and call the simulation methods to actually simulate and record the model.  The Matlab file script.m performs all of these actions in one script.  Most of the parameter values used are stored in this script and copied to the corresponding fields in the Network object.


1 - Quick-start guide for running the model


If not already done, unpack the contents of the .zip file containing these instructions to a common folder, and in Matlab navigate to this folder.  Run script.m by typing 'script' into the interpreter.  This will

a) create a directory 'runResults' where the recorded simulation data is to be extracted;

b) load the precompiled Java file SimJNSVer14-6.jar into the Matlab classpath;

c) instantiate an object of type Network;

d) copy parameter values into the Network object;

e) call model construction routines of the Network class;

f) make small modifications to connectivity not included in Network class;

g) call the script runAndRecordJavaSim, which calls the Network run methods, periodically records Network state, and finally saves the result in a Matlab .mat file in the folder runResults.

The resulting run[x].mat file contains the variables spikesExc and spikeTimesExc, which are 1 x n vectors containing the excitatory neuron indexes and times corresponding to each excitatory spike that occurred in the simulation recording window.  The variable numSpikesExc contains the number of excitatory spikes that occured, so executing the lines

figure;
scatter(spikeTimesExc(1:numSpikesExc),spikesExc(1:numSpikesExc),'r.')

will produce a raster plot of all excitatory spikes during simulation.  The file run[x].mat also contains corresponding vectors for inhibitory spikes.


2 - Compiling Java source code


While the precompiled Java source code is already supplied in the .zip file as SimJNSVer14-6.jar, the uncompiled source code is also included in the src folder.  To reproduce the precompiled java code, compile the java source code in the src folder into the jar file SimJNSVer14-6.jar, and include the dependencies on the external jar files lib/cern.jar and lib/jmatio.jar.
