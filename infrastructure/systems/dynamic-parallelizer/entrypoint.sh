#!/bin/bash
base=$(dirname $0)
source ${base}/.venv/bin/activate
exec "$@"
