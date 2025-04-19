#!/usr/bin/env python3

from pathlib import Path
from time import sleep

# we use python here because python bytes does not have asymptotic overhead (almost exactly 1G)

path = Path(__file__).resolve().parent.parent / 'inputs' / '1G.zeros'
one_gig = path.read_bytes()
sleep(1)
