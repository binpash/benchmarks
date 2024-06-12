#!/bin/bash

set -e

PYTHON=${PYTHON:-`which python`}
DIR=$PWD
TMP=$DIR/tmp
SCRIPTS=$DIR/scripts

# Ideally, we'll move on to piping rather than writing to a file
MODEL=$TMP/model.obj
X=$TMP/X_train.obj
y=$TMP/y_train.obj
CLASSES=$TMP/classes.obj
DUAL=false # should be converted to bool inside script
MAX_SQ_SUM=$TMP/max_squared_sum.obj
WARM_COEF=$TMP/warm_start_coef.obj
C_=$TMP/C_.obj

# TODO: Try this out on a larger dataset
# TODO: Benchmark each phase

# Generating model & samples
$PYTHON $SCRIPTS/gen_model.py 100
$PYTHON $SCRIPTS/gen_samples.py

# Validity checking functions
# These functions just check to make sure that the input is valid. 
# If not they will raise an error. Otherwise, they do not mutate the data.
$PYTHON $SCRIPTS/check_solver.py $MODEL
penalty=$($PYTHON $SCRIPTS/penalty.py $MODEL)
$PYTHON $SCRIPTS/val_data.py $MODEL $X $y 
$PYTHON $SCRIPTS/classes.py $MODEL $y # This should return a classes with just the unique classes in y
multiclass=$($PYTHON $SCRIPTS/check_multiclass.py $MODEL)

# TODO: Benchmark each step of the pipeline
# Make a modified pipeline where each step writes its output to a file

# Calculations functions
$PYTHON $SCRIPTS/rownorm.py $X
n_classes=$($PYTHON $SCRIPTS/reshape_classes.py $MODEL $CLASSES)
$PYTHON $SCRIPTS/warm_start.py $MODEL $multiclass $n_classes # pipes coefficients

# RSA1 dataset has 104 classes
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 1
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 2
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 3
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 4
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 5
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 6
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 7
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 8
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 9
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 10
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 11
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 12
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 13
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 14
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 15
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 16
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 17
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 18
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 19
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 20
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 21
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 22
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 23
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 24
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 25
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 26
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 27
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 28
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 29
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 30
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 31
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 32
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 33
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 34
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 35
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 36
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 37
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 38
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 39
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 40
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 41
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 42
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 43
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 44
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 45
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 46
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 47
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 48
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 49
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 50
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 51
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 52
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 53
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 54
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 55
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 56
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 57
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 58
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 59
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 60
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 61
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 62
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 63
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 64
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 65
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 66
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 67
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 68
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 69
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 70
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 71
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 72
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 73
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 74
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 75
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 76
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 77
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 78
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 79
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 80
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 81
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 82
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 83
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 84
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 85
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 86
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 87
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 88
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 89
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 90
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 91
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 92
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 93
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 94
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 95
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 96
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 97
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 98
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 99
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 100
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 101
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 102
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 103
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 104

$PYTHON $SCRIPTS/zip_coef.py $MODEL $n_classes
$PYTHON $SCRIPTS/adjust_coef.py $MODEL $X $multiclass $n_classes result/trained_model.obj
