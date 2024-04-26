#!/bin/bash

scripts=(
    1_1.sh
    2_1.sh
    2_2.sh
    3_1.sh
    3_2.sh
    3_3.sh
    4_3.sh
    4_3b.sh
    6_1.sh
    6_1_1.sh
    6_1_2.sh
    6_2.sh
    6_3.sh
    6_4.sh
    6_5.sh
    6_7.sh
    7_1.sh
    7_2.sh
    8.2_1.sh
    8.2_2.sh
    8.3_2.sh
    8.3_3.sh
    8_1.sh
)


for script in ${scripts[@]}
    do
        ./$script
done
