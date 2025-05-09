echo hi
for i in 1 2; do
    echo hi1
    for j in 1 2 3; do
        echo hi2
        sleep 0.2
        echo hi3
    done
    echo hi4
done
echo hi5

echo hi6
for i in 1 2 3; do
    echo hi7
    for j in 1 2; do
        echo hi8
        sleep 0.2
        echo hi9
    done
    echo hi10
done
echo hi11


## Future loop tests must include:
## 1. A single loop with a single command without anything else 
## 2. Multiple commands in the same loop
## 3. Nested loops
## 4. Commands before and after a loop
