#!/usr/bin/env python3

import datetime
import io
import matplotlib.pyplot as plt
import numpy as np
import sys

# Read data from stdin
data = sys.stdin
year = sys.argv[1]

# Transform data from string to datetime and float
to_dt = lambda x: datetime.datetime.strptime(f"{year}-{x}", "%Y-%m-%d")
data = [ line.strip().replace("\r", "") for line in data ]
data = [ (to_dt(date), float(_min), float(_max), float(temp)) for date, _min, _max, temp in [ line.split() for line in data if line ] ]
data = { key: [ x for x in values ] for key, values in zip(["date", "min", "max", "temp"], zip(*data)) }

# List of dates
dates = data["date"]
# Max temperature across all years
historic_max = data["max"]
# Min temperature across all years
historic_min = data["min"]
# Temperature for the year
temps = data["temp"]

# Plotting
plt.figure(figsize=(10, 5))

# Labels and titles
plt.xlabel(f"Days in {year}")
plt.ylabel("Temperature (°F)")
plt.title(f"Temperature Records for {year}")

# Save the plot so far
plt.savefig(f"{year}_onlyaxes.png")

# Plot the mimimum and maximum temperature range
plt.fill_between(dates, historic_min, historic_max, color="tan", alpha=0.3, label='Normal range')

# Save the plot so far
plt.savefig(f"{year}_range.png")

# Plot given year's temperature
plt.plot(dates, temps, color="black", label=f"{year} temperature")

# Save the plot so far
plt.savefig(f"{year}_temperature.png")

# Legend
plt.legend(loc="upper right")

# Show the plot
plt.grid(True)

# Save the plot so far
plt.savefig(f"{year}_final.png")

buf = io.BytesIO()
plt.savefig(buf, format='png')
plt.close()
sys.stdout.buffer.write(buf.getvalue())
