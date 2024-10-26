# Ben.sh

> _It's a suite that has benchmarks in it._

## Docker

```sh
# Build the container
$ docker build -t bensh .

# Run the container
$ docker run -it bensh

# For development, mount the benchmarks directory
docker run -it -v "$(pwd):/benchmarks" bensh
```
