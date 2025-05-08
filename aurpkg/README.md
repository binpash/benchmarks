## aurpkg

This benchmark converts a list of package names into `pacscript` files by fetching and parsing their corresponding `PKGBUILD` files.

### Inputs

- `inputs/packages`: A newline-delimited list of package names.

### Running

For each package name in the input file:

1. The corresponding `PKGBUILD` file is retrieved.
2. The following fields are extracted:
   - `pkgname`
   - `pkgver`
   - `pkgdesc`
   - `url`
   - `source`
   - `sha256sums`
   - `depends`
   - `makedepends`
3. These values are written into a `pacscript` containing stub `build()` and `package()` functions.
4. The result is saved to `outputs/<package>/<package>.pacscript`.

### Validation

Correctness is determined by computing the SHA-256 hash of each generated pacscript file and comparing it against a reference hash stored in hashes/.
