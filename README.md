# bitvector_hw2

The github link for this project is: https://github.com/bgoldsc1/bitvector_hw2

The project is written in C++ and each task is completed in the file main.cpp. This has the following classes: Bitvector, Rank, Select, and Sparse. Note that, in order to create a Bitvector, one must pass a string representation of the Bitvector to the constructor for the Bitvector class. For example, to create the Bitvector 01010110 and store it in the variable b, one would write

Bitvector b = Bitvector("01010110");

Note that, in actuality, the Bitvector class at time of construction appends a certain number of 0's at the end of the given string so that the length of a Bitvector object is always an even power of 2, i.e. can be written in the form 2^(2N) where N is a positive integer. This increases the size of the bitvector by at most a constant factor (by at most 4 times, in fact). 

Note that I use the compact vector third party library in my implementation of the Bitvectors. I had to use the version in Professor Patro's pufferfish github with certain sections commented out. The versions in my github repository are the versions that work with my implementation for this project. 

Task 1

The Rank class fulfills the requirements of the first task. An object in this class can be created from a Bitvector object by passing a pointer to a Bitvector object into the constructor for the Rank class. For example, one would write

Rank r = Rank(&b);

The remaining functions are exactly as specified in the spec, although the filename is passed to the save and load functions by value rather than by reference.

Task 2

The Select class fulfills the requirements of the second task. An object in this class can be created from a Rank object by passing a pointer to a Rank object into the constructor for the Select class. For example, one would write

Select s = Select(&r);

The remaining functions are exactly as specified in the spec, although the filename is passed to the save and load functions by value rather than by reference.

Task 3

The Sparse class fulfills the requirements of the third task. An object in this class can be created from a size, which constructs a Bitvector of 0's of the specified size. For example, one would write

Sparse sp = Sparse(10);

in order to construct an empty sparse array of size 10. Note that the function which is called create in the spec here is simply implemented using this form of the constructor for the Sparse class. 

The remaining functions are exactly as specified in the spec, although the filename is passed to the save and load functions by value rather than by reference.
