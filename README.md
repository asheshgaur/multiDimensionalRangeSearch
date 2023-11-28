# multiDimensionalRangeSearch
Program for creating multidimensional range search tree.

Run make to install the program with the dependencies mentioned in the Makefile.
The program should be run like this:

./main {Rules_File} {Max_Rules_Per_Node} {Brute_Classification - 0/1}

Max_Rules_Per_Node specifies the maximum number of rules to be stored in a single node of the tree. It is similar to 'binth' defined in HiCuts. An integer value will be provided.

Brute_Classification is optional and was used for correctness analysis for the project. Use 0 to bypass it.