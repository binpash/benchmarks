#!/bin/bash

# Generate a large number by multiplying prime numbers.
# This is just a simplistic way to generate a large number.
# Note: This is not intended to be an efficient prime number generator.

# Generate a large number (e.g., product of first few primes raised to a power)
large_number=2
filename=$1

touch "$filename"

for i in {3..9}; do
  large_number=$((large_number * $i))
done
echo $large_number
# Save this large number into a file
echo $large_number > "$filename"

