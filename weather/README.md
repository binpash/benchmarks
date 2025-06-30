## weather

The `max-temp` part of the benchmark processes NOAA weather data to compute basic temperature statistics across multiple years.
The `tuft-weather` part re-creates Edward Tuft's famous weather plot using the same, or similar, data.

### Inputs

- `inputs/temperatures.txt`: A plain-text file containing fixed-width weather records.

### Running

The input file is processed to compute summary statistics for maximum, minimum, and average temperature across all data points.
Each value is computed by extracting a fixed character range from each line, filtering out sentinel values, and applying standard Unix text processing tools.

### Validation

Correctness is determined by comparing each output file to its corresponding reference file in `correct-results/`.
