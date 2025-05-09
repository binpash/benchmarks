# a=$(free -m | awk '/Mem:/ { print $4 }')
# b=$(python3 -c "print(${a}*0.7)")

# base=${HOME}/tmp
# ${base}/alloc_echo 3 2 hello1
# ${base}/alloc_echo 3 2 hello2
# ${base}/alloc_echo 1 4 hello3
# ${base}/alloc_echo 3 4 hello4
# ${base}/alloc_echo 1 4 hello5
# ${base}/alloc_echo 3 4 hello6

# tail --byte 1073741824 /dev/zero > x1
# tail --byte 1073741824 /dev/zero > x2
# tail --byte 1073741824 /dev/zero > x3
# tail --byte 1073741824 /dev/zero > x4
# tail --byte 1073741824 /dev/zero > x5
# tail --byte 1073741824 /dev/zero > x6
# tail --byte 1073741824 /dev/zero > x7
# tail --byte 1073741824 /dev/zero > x8
# tail --byte 1073741824 /dev/zero > x9

# mv x y
# mv y z
# mv z w
# mv w u
# mv u t
# mv t x

a=$(echo 10)
if [ $a -eq 11 ]; then
    echo oh
    exit
else
    echo ho
    exit
fi
