#!/bin/bash

PYTHON="python3"
OUT=${OUT:-$PWD/result}
TMP=${TMP:-$PWD/tmp}
#export tmp to env
export TMP
SCRIPTS=${SCRIPTS:-$PWD/scripts}

# Ideally, we'll move on to piping rather than writing to a file
MODEL=$TMP/model.obj
X=$TMP/X_train.obj
y=$TMP/y_train.obj
CLASSES=$TMP/classes.obj
DUAL=false # should be converted to bool inside script
MAX_SQ_SUM=$TMP/max_squared_sum.obj
WARM_COEF=$TMP/warm_start_coef.obj
C_=$TMP/C_.obj

echo $PYTHON >&2
echo "DIR: $DIR" >&2
echo "SCRIPTS: $SCRIPTS" >&2
echo "MODEL: $MODEL" >&2
echo "X: $X" >&2
echo "y: $y" >&2
echo "CLASSES: $CLASSES" >&2
echo "DUAL: $DUAL" >&2
echo "MAX_SQ_SUM: $MAX_SQ_SUM" >&2
echo "WARM_COEF: $WARM_COEF" >&2
echo "C_: $C_" >&2

# Generating model & samples
$PYTHON $SCRIPTS/gen_model.py 100
$PYTHON $SCRIPTS/gen_samples.py "$@"

# Validity checking functions
# These functions just check to make sure that the input is valid. 
# If not they will raise an error. Otherwise, they do not mutate the data.
$PYTHON $SCRIPTS/check_solver.py $MODEL
penalty=$($PYTHON $SCRIPTS/penalty.py $MODEL)
$PYTHON $SCRIPTS/val_data.py $MODEL $X $y 
$PYTHON $SCRIPTS/classes.py $MODEL $y # This should return a classes with just the unique classes in y
echo "$PYTHON $SCRIPTS/check_multiclass.py $MODEL" >&2
multiclass=$($PYTHON $SCRIPTS/check_multiclass.py $MODEL)
echo "------" >&2
# Make a modified pipeline where each step writes its output to a file

# Calculations functions
$PYTHON $SCRIPTS/rownorm.py $X
n_classes=$($PYTHON $SCRIPTS/reshape_classes.py $MODEL $CLASSES)
$PYTHON $SCRIPTS/warm_start.py $MODEL $multiclass $n_classes # pipes coefficients

# Covtype dataset has 7 classes
echo "WARM_COEF: $WARM_COEF" >&2
echo "MAX_SQ_SUM: $MAX_SQ_SUM" >&2

echo "multiclass: $multiclass" >&2
echo "penalty: $penalty" >&2
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 1
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 2
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 3
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 4
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 5
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 6
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 7

$PYTHON $SCRIPTS/zip_coef.py $MODEL
$PYTHON $SCRIPTS/adjust_coef.py $MODEL $X $multiclass $n_classes $RESULT/trained_model.obj
