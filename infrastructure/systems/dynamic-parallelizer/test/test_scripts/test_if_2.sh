a=hello
if [ -z $a ]; then
    echo what
fi

if [ ! -z $a ]; then
    echo who
fi
