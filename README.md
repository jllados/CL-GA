Genetic Algorithm capable of obtaining the best possible combinations of constraints for a consistency library
CL-GA
=========
GAlib is a set of C++ genetic algorithm objects.
The library includes tools for using genetic algorithms to do optimization in any C++ program.


Prerequisites
--------------
GAlib compilation requires the following tools installed on your system ``make``, ``gcc-c++`` and ``t_coffee``. 


Compile 
--------
Download the git repository on your computer.
    
Make sure you have installed the required dependencies listed above. 
When done, move in the project root folder and enter the following commands:     
    
    $ cd examples
    $ make
    

The binary will be automatically generated in the folder.


How to use
--------

    ./GA_TC_lib output_name length_of_the_gen


Output files
--------

    _bog.dat (best gen for each flushFrequency)
    .gen (indexes of the gen)
    .lib (library in T-Coffee format of the .gen)
    .msf (alignment output)
