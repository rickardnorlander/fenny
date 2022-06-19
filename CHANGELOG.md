# 0.2.0

* Added binary search for lexicographic version of 2d and 3d.
* Renamed lexicographic functions from so\_x to x\_so.
* Renamed lexicographic functions from x\_linear to x\_lex.
* Changed the output of `get_root_p1`, for inputs of the form (2^n) - 1. Note that the binary search routines worked correctly with the previous value as well.
* Changed z dimension to be least significant (this is an internal change).
* Reduced cpu and memory usage for the lexicographic functions.
* Added more testing.
* Added rustdocs.
