dynamic_exit ()
{
    echo about to exit
    exit 5
    echo not reached
}

echo before call
dynamic_exit
echo not reached
