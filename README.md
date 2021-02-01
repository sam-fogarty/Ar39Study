# Ar39Study

This is main script and some supporting scripts used for the Argon-39 beta-decay analysis project.
This uses basic pattern recognition to find point-like beta-decay candidates. 
The project utilizes Art/ROOT, Gallery, and LArSoft, and works on ProtoDUNE and MicroBooNE. 

**Set up LArSoft and Gallery** (setup scripts provided) and build the tool with

```$ make```


## Ar39Study Tool

The Ar39Study tool requires the *param.txt* and *filelist.txt* files be present. A *param.txt*
is provided, a *filelist.txt* must be generated in the same path as *param.txt*. The filelist is 
made of the ***full paths and names*** of the data to be analyzed. 

To run Ar39Study with parameters as set in *param.txt*, 

```$ ./Ar39Study```

Number of event readouts analyzed can be changed in *param.txt*. A specific range of readouts can 
be selected using the following syntax

```$ ./Ar39Study [first event] [last event]```

An arbitrary file can be fed to Ar39Study like this

```$ ./Ar39Study /path/to/arbitrary/file/filename```

If given a file this way, Ar39Study will use the maximum events from *param.txt*. This is useful 
for scripts used to run this tool on the grid to analyze many files in parallel. 


## MicroBooNE Specific Features

### Point-Like 3D Position

There is an experimental feature that uses MicroBooNE geometry to try to infer 3D position of point-like 
signals using the induction wires as well as the collection wires. This requires the MicroBooNE 
*ChannelWireGeometry_V2.txt* file with the top text ***removed***. It should only be a table of numbers. 
Then the second required file, CollectionIntersection.txt can be generated with /scripts/MicroBoonE/GeomScript. 

This feature is expressed in the Confirm_Candidates() function in *PointSignalFunctions.h*.

