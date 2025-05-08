## covid-mts

This benchmark analyzes CSV logs of vehicle movement and produces multiple statistical summaries based on time, vehicle, and date fields.

### Inputs

- `inputs/in.csv`: A CSV file containing timestamped records of vehicle activity.

### Running

For each of five analysis scripts:

1. The input file is processed to extract fields such as date, time, vehicle ID, and route.
2. Specific aggregations are performed:
   - `1.sh: Counts the number of vehicles active per day.
   - `2.sh: Counts the number of days each vehicle was active.
   - `3.sh: Counts the number of hours each vehicle was active.
   - `4.sh: Counts the number of distinct hours recorded per day.
   - `5.sh: Tabulates the number of active hours per vehicle across all dates.
3. The output of each script is saved to `outputs/outputs/<script_number>.out`.

### Validation

Correctness is determined by computing the MD5 hash of each output file and comparing it against a reference hash stored in `hashes/`.

### Refenreces

- https://insidestory.gr/article/noymera-leoforeia-athinas
