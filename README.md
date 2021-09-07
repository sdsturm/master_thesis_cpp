# Master's Thesis C++ Codes

## Structure

- `./lib`: Main library (`mthesis-lib`)
- `./lib/test`: Unit tests for main library
- `./production`: Executables to produce results for the actual thesis
- `./examples`: Example applications to demonstrate the use of `mthesis-lib`
- `./submodules`: External modules which are typically not available by the 
system package manager.

## Dependencies

See `CMakeLists.txt` files at top level `./` and in `./lib` and `./submodules`.

Note: [gnuplot-iostream](https://github.com/dstahlke/gnuplot-iostream.git)
needs the [gnuplot](http://www.gnuplot.info/) application to be found by the
system shell.
