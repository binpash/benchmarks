#!/bin/bash

seq 1 "$1" | xargs -n1 -I{} "$TRY" -y echo x
