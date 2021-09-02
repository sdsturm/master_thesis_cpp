# Propagating Plane-Wave Representations of Half-Space Green's Functions

## Structure

- `./lib`: Main library (`mthesis-lib`)
- `./production`: Executables to produce results for the actual thesis
- `./results`: Generated data output files (ignored by Git)
- `./examples`: Example applications to demonstrate the use of `mthesis-lib`

## Dependencies

### System

- [Boost](https://www.boost.org/) (only the header-only parts and `Boost::test`
for the main library)
- [Armadillo](http://arma.sourceforge.net/)
- openBLAS (to link against instead of armadillo runtime library)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)

### Submodules Dependencies

#### [complex_bessel](https://github.com/joeydumont/complex_bessel.git)

- Fortran compiler
- HDF 5
- Google Test (optional to run tests)

#### [gnuplot-iostream](https://github.com/dstahlke/gnuplot-iostream.git)

- Installation of [gnuplot](http://www.gnuplot.info/) at default location
- Boost binary components `Boost::iostreams`, `Boost::system`,
`Boost::filesystem`
