# Koala

> _It's a suite that has benchmarks in it._

## Benchmarks

| Benchmark    | Description                                             |
| ---------    | -----------                                             |
| aurpkg       |                                                         |
| bio          | Bioinformatics.                                         |
| covid-mts    | COVID-19 multivariate time series.                      |
| file-enc     | File encoding.                                          |
| log-analysis | Log analysis.                                           |
| max-temp     | Maximum temperature.                                    |
| media-conv   | Media conversion.                                       |
| nlp          | Natural language processing.                            |
| oneliners    | One-liners.                                             |
| riker        |                                                         |
| sklearn      | Machine learning.                                       |
| teraseq      |                                                         |
| tuft-weather | Tuft weather.                                           |
| uniq-ips     | Unique IPs.                                             |
| unix50       | Unix 50.                                                |
| web-index    | Web index.                                              |

## Docker Instructions

```sh
# Build the container
$ docker build -t koala .

# Run the container
$ docker run --cap-add NET_ADMIN --cap-add NET_RAW -it koala

# For development, mount the benchmarks directory
docker run --cap-add NET_ADMIN --cap-add NET_RAW -it -v "$(pwd):/benchmarks" koala
```


⠀⠀⠀⣠⠖⠚⠉⠙⠲⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠾⠋⠉⠑⠢⣄⠀⠀⠀⠀⠀
⠀⢀⡞⠁⠀⠀⠀⠀⠀⠈⠓⠦⠤⠤⠴⠖⠒⠉⠉⠉⠉⠉⠉⠉⠒⠲⠦⠤⠤⠴⠟⠁⠀⠀⠀⠀⠀⠈⢷⡀⠀⠀⠀
⠀⣾⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣧⠀⠀⠀
⢸⡇⢀⡴⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢦⠀⢹⣆⠀⠀
⠈⢧⡞⣴⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⣶⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠠⣤⢧⡾⠙⠀⠀
⠀⠈⢿⡇⠀⠀⠀⢠⠇⠀⠀⠀⠀⠀⠀⠀⠀⣾⣣⣿⣿⣿⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⠀⠀⠀⠀⢹⡾⠀⠀⠀⠀
⠀⠀⠘⣇⠀⠀⢀⡾⠀⠀⠀⠀⣾⣿⠀⠀⢸⣯⣿⣿⣿⣿⡇⠀⠀⣿⣷⠄⠀⠀⠀⠀⠈⣷⠀⠀⠀⣼⠁⠀⠀⠀⠀
⠀⠀⠀⠈⠳⠤⣼⠁⠀⠀⠀⠀⠈⠁⠀⠀⠸⣿⣿⣿⣿⣿⡇⠀⠀⠈⠉⠀⠀⠀⠀⠀⠀⠘⣧⠤⠞⠁⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠻⠿⠿⠿⠟⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣹⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢿⡄⠀⠀⠀⠀⠀⠀⠀⠀⠘⠲⠴⠟⠦⠴⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⡟⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠹⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⢦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣶⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⠋⠁⠀⠀⠀⠀⠀⢀⡖⠚⣻⠀⠀⠀⠀⠀⠀⠀⠈⢻⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡏⠀⠀⠀⣾⠉⠙⢲⣼⡤⠞⠉⠀⠀⠀⢐⡆⠀⠀⠀⠈⣯⢧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠇⠀⠀⠀⠈⣩⠉⠉⠙⣆⠀⠀⠀⠀⠀⣸⠇⠀⠀⠀⠀⢹⡈⢷⡄⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⠀⠀⠀⠀⢰⡏⠀⠀⠀⠸⡆⠀⠀⣠⡶⠋⠀⠀⠀⠀⠀⢸⡇⠀⢿⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⡏⠀⠀⠀⠀⢸⡇⠀⠀⠀⢀⡿⠶⠋⠉⠀⠀⠀⠀⠀⠀⠀⣼⠃⠀⠘⣗⠲⣆⠀⠀⠀⠀⠀
⠀⠀⠀⠀⢠⡖⠋⠙⠲⡇⠀⠀⠀⠀⢸⠀⣰⠤⢥⣿⠤⠤⠄⠀⠀⠀⠀⠀⣀⣤⠞⠁⠀⠀⠀⢹⡆⠛⢶⡄⠀⠀⠀
⠀⠀⠀⠀⡟⠀⠀⠀⠀⡇⠀⠀⠀⠀⢸⣸⠃⠀⠀⠙⢷⣖⣃⣠⡤⠴⠖⠋⠉⠀⠀⠀⠀⠀⠀⠈⣷⠀⡼⠃⠀⠀⠀
⠀⠀⠀⢸⡇⠀⠀⠀⠀⡇⠀⠀⠀⠀⢸⡟⠀⠀⠀⠀⠈⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⣠⠏⠀⠀⠀⠀
⠀⠀⠀⠈⣧⠀⠀⠀⠀⣷⠀⠀⠀⠀⣼⣇⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡿⠋⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠘⣧⠀⠀⠀⢿⣀⡀⠀⠀⡿⢻⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⡴⠟⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠈⠳⠤⠤⠼⢷⣇⣸⡾⠓⠚⠳⣄⣀⠀⠀⠀⠀⣀⣀⣀⣠⡤⠤⠴⠖⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠉⠉⠉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
