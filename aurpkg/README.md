## aurpkg

This benchmark reads a list of package names, fetches the AUR package build scripts (`PKGBUILD`) and builds them using the `makedeb` tool.

### Inputs

- `inputs/packages`: A newline-delimited list of package names.

### Running

For each package name in the input file:

1. The corresponding `PKGBUILD` file is retrieved.
2. The script attempts to import any GPG keys declared in `validpgpkeys=()` within the `PKGBUILD` (if present).
3. The `makedeb -d` command is run to simulate the creation of a Debian package.
4. The results are stored as:
- `outputs/<package>/PKGBUILD` – the downloaded PKGBUILD
- `outputs/<package>.txt` – stdout and stderr log from `makedeb`

If the script is executed as root, it uses `runuser` to execute `makedeb` as a non-root user (`user`).

### Validation

Correctness is determined by checking whether the log file contains the string:
```bash
Finished making
```
This is performed using:
```bash
grep "Finished making" outputs/<package>.txt
```
If the string is found for all packages, the run is considered successful.


