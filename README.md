# dglap90

**dglap90** is a package of Fortran 90 code for peforming DGLAP evolution
on unpolarized parton distributions.
Currently, it is implemented only at leading order.
It uses a brute force method, similar to that employed by Kumano and Miyama.

**dglap90** reads and writes CSV (comma separated values) files,
with each column identified by a key in the header.
In this way, the order of the columns for the parton distributions
does not matter.
The code will also not fail if a column (other than the x column) is missing,
as the missing parton will be assigned zero values.

## Building the test program

```bash
$ git clone https://github.com/adamfreese/dglap90.git /path/to/source
$ cd /path/to/build
$ cmake /path/to/source
$ make
```

## Using dglap90 in other code

The file dglap.f90 contains the DGLAP module, but depends on modules
contained in basics.f90, dataframes.f90, and qcd.f90.
One could compile all of these files together to get a single library.
If compiled separately, when linking, note that everything else depends
on basics, and dglap depends on everything else.

In any code that is to use the DGLAP module, one should have:
```fortran90
use dglap
```

There are two ways to initialize a PDF grid to be evolved.
One is to call
```fortran90
call dglap_init( opt_Nx=Nx, opt_xmin=xmin, opt_xmax=xmax, opt_Q2step=Q2step )
call dglap_set_pdf( pdf_function, Q2_value )
```
where `pdf_function` takes a double-precision real as input and gives
an array of double-precision reals, with indices from -6 to 6, as output.
Note that up is associated with the index 1, and down with the index 2.
The call to `dglap_init` is optional,
and all of its arguments are also optional.

The second method is to read a CSV file using
```fortran90
call dglap_read_csv( '/path/to/csv/file', Q2_value )
```
The CSV file should have at least a column for *x*,
with the label `x` in the header.
The parton labels are
`u`, `d`, `s`, `c`, `b`, `t`, `g`,
`ub`, `db`, `sb`, `cb`, `bb`, and `tb`.
The order of the columns in the CSV file does not matter.
Any partons missing from the CSV file will be assigned 0 values.

Both initialization methods require giving a Q-squared value (in GeV-squared).

To evolve a PDF grid, call
```fortran90
call dglap_evolve( new_Q2_value )
```

To save the current PDF grid to a file, call
```fortran90
call dglap_write_csv( '/path/to/csv/file' )
```

**dglap90** dynamically allocates memory to be used internally.
To explicitly deallocate the memory when done, call
```fortran90
call dglap_cleanup()
```

## References

The test program uses a leading order PDF from the CJ15 set.
More information can be found
[on their website](https://www.jlab.org/theory/cj/pdfs.html).
See also:

- A. Accardi, L. T. Brady, W. Melnitchouk, J. F. Owens, and N. Sato,
[Phys. Rev. **D93** (2016), 114017](http://inspirehep.net/record/1420566)

The use of brute force integration to solve the DGLAP equation is very
similar to that used by Kumano and Miyama, *cf.*,

- M. Miyama and S. Kumano,
[Comput. Phys. Commun. 94 (1996) 185-215](http://inspirehep.net/record/397989)

## License

**dglap90** is available under the GNU GPLv3.
See the file COPYING for the license.
