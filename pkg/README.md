## pkg

This benchmark reads a list of package names, fetches the AUR package build scripts (`PKGBUILD`) and builds them using the `makedeb` tool.
It also performs inference of security permissions on a number of npm packages.

### Inputs

- `inputs/packages`: A newline-delimited list of package names.

### Validation

For the `pacaur` part of the benchmark, correctness is determined by checking whether the log file contains the string "Finished making".
If the string is found for all packages, the run is considered successful.

For the permission inference part, correctness is determined by that the inferred permissions match the pre-computed ones through their hash.
