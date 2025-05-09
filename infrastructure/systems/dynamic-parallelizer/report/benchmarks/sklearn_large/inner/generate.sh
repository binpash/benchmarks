#!/bin/bash

# Script to generate lines that launch individual regressors. We need n regressors for n class classifications
n_classes=104
for i in `seq 1 $n_classes`
do
	echo "\$PYTHON \$SCRIPTS/parallel.py \$MODEL \$X \$y \$C_ \$WARM_COEF \$MAX_SQ_SUM \$multiclass \$penalty $i"
done
