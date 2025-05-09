#!/bin/bash

export ORCH_TOP=${ORCH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
export WORKING_DIR="$ORCH_TOP/test"
export TEST_SCRIPT_DIR="$WORKING_DIR/test_scripts"
export MISC_SCRIPT_DIR="$WORKING_DIR/misc"

echo "==================| Scheduler Tests |==================="
echo "Test directory:               $WORKING_DIR"
echo "Test script directory:        $TEST_SCRIPT_DIR"

## Set the DEBUG env variable to see detailed output
DEBUG=${DEBUG:-0}

bash="bash"
## Debug needs to be set to 2 because otherwise repetitions cannot be checked
orch="$ORCH_TOP/pash-spec.sh -d 2"
# Generated test scripts are saved here
test_dir_orch="$ORCH_TOP/test/test_scripts_orch"
test_dir_bash="$ORCH_TOP/test/test_scripts_bash"
# Test script output is saved here
output_dir_orch="$ORCH_TOP/test/output_orch"
output_dir_bash="$ORCH_TOP/test/output_bash"
# Results saved here
output_dir="$ORCH_TOP/test/results"

echo "Bash scripts saved at:       $test_dir_bash"
echo "Orch scripts saved at:       $test_dir_orch"
echo "Results saved at:            $output_dir"
echo "========================================================"

# Clear previous test results
rm -rf "$output_dir"
mkdir -p "$output_dir"
touch "$output_dir/result_status"

cleanup()
{
    # clear Riker's cache
    rm -rf ./.rkr
    rm -rf "$output_dir_orch"
    rm -rf "$output_dir_bash"
    mkdir "$output_dir_orch"
    mkdir "$output_dir_bash"
}

test_repetitions()
{
    local repetitions="$1"
    local exec_log_file="$2"
    if [ "${#repetitions}" -eq 1 ]; then
        result=`python3 $WORKING_DIR/parse_cmd_repetitions.py "--total" $exec_log_file`
        if [ "$repetitions" != "$result" ]; then
            echo " (!) Reps (total) not optimal: Expected: $repetitions | Got: $result" 1>&2
            return 1
        fi
    else
        result=`python3 $WORKING_DIR/parse_cmd_repetitions.py "--detailed" $exec_log_file`
        if [ "$repetitions" != "$result" ]; then
            echo " (!) Reps (detailed) not optimal: Expected: $repetitions | Got: $result" 1>&2
            return 1
        fi
    fi
}

run_test()
{
    cleanup
    local test=$1
    shift # Move past the test name to process potential execution time limit and repetitions
    local execution_time_limit=$1
    shift # Move past the execution time limit (could be empty if not specified)
    local repetitions=("$@") # Remaining arguments are treated as repetitions

    if [ "$(type -t $test)" != "function" ]; then
        echo "$test is not a function!   FAIL"
        return 1
    fi

    printf "Running %-35s" "$test..."

    # Bash execution
    failure_reason="" # Variable to store the reason of failure
    export test_output_dir="$WORKING_DIR/output_bash"
    $test "$bash" "$TEST_SCRIPT_DIR" "$test_output_dir" "$execution_time_limit" "${repetitions[@]}" > "$test_output_dir/stdout" 2> /dev/null
    test_bash_ec=$?

    # Orch execution
    export test_output_dir="$WORKING_DIR/output_orch"
    TIMEFORMAT='%3R'; time ( $test "$orch" "$TEST_SCRIPT_DIR" "$test_output_dir" "$execution_time_limit" "${repetitions[@]}" > "$test_output_dir/stdout" 2> "$test_output_dir/test_stderr" ) 2> "$test_output_dir/time_output"
    test_orch_ec=$?

    # If debug is set, print the stderr
    if [ "$DEBUG" -gt 1 ]; then
        cat "$test_output_dir/test_stderr"
    fi

    local actual_execution_time=$(cat "$test_output_dir/time_output")

    # Check diffs
    diff --recursive --exclude="time_output" --exclude="test_stderr" "$WORKING_DIR/output_bash/" "$WORKING_DIR/output_orch/"
    test_diff_ec=$?

    test_execution_time_ec=0
    # Execution time testing
    if [ ! -z "$execution_time_limit" ] && (( $(echo "$actual_execution_time > $execution_time_limit" | bc -l) )); then
        failure_reason="(?) Execution time exceeded: Limit: $execution_time_limit sec, Actual: $actual_execution_time sec"
        test_execution_time_ec=1
        echo "$test $failure_reason" >> $output_dir/result_status
    fi

    # Test repetitions if specified
    if [ "${#repetitions[@]}" -gt 0 ]; then
        test_repetitions "${repetitions[*]}" "$stderr_file"
        test_repetitions_ec=$?
        if [ $test_repetitions_ec -ne 0 ]; then
            failure_reason=${failure_reason:-"(?) Repetitions mismatch"}
        fi
    fi

    # Evaluate results
    if [ $test_diff_ec -ne 0 ]; then
        failure_reason=${failure_reason:-"(!) Output mismatch"}
    fi

    if [ $test_bash_ec -ne $test_orch_ec ]; then
        failure_reason=${failure_reason:-"EC mismatch [$test_bash_ec-$test_orch_ec]"}
    fi

    if [ -n "$failure_reason" ]; then
        echo -en "FAIL | Time: $actual_execution_time sec | $failure_reason\n"
    else
        echo -en "OK   | Time: $actual_execution_time sec\n"
    fi

    if [ "${test_diff_ec:-0}" -ne 0 ] || [ "${test_repetitions_ec:-0}" -ne 0 ] || [ "${test_execution_time_ec:-0}" -ne 0 ]; then
        echo "$test: FAIL" >> $output_dir/result_status
    else
        echo "$test: OK" >> $output_dir/result_status
    fi
}

test_single_command()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > "$3/in1"
    $shell $2/test_single_command.sh
}

test1_1()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > "$3/in1"
    $shell $2/test1_1.sh
}

test1_2()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > "$3/in1"
    $shell $2/test1_2.sh
}

test1_3()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > "$3/in1"
    $shell $2/test1_3.sh
}

test2_1()
{
    local shell=$1
    $shell $2/test2_1.sh
}

test2_2()
{
    local shell=$1
    $shell $2/test2_2.sh
}

test2_3()
{
    local shell=$1
    $shell $2/test2_3.sh
}

test3_1()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in2
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in3
    $shell "$2/test3_1.sh"
}

test3_2()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in2
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in3
    $shell "$2/test3_2.sh"
}

test3_3()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in2
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in3
    $shell "$2/test3_3.sh"
}

test4_1()
{
    local shell=$1
    echo 'hello1' > "$3/in1"
    echo 'hello2' > "$3/in2"
    $shell "$2/test4_1.sh"
}

test4_2()
{
    local shell=$1
    echo 'hello1' > "$3/in1"
    echo 'hello2' > "$3/in2"
    $shell "$2/test4_2.sh"
}

test4_3()
{
    local shell=$1
    echo 'hello1' > "$3/in1"
    echo 'hello2' > "$3/in2"
    $shell "$2/test4_3.sh"
}

test5_1()
{
    local shell=$1
    echo 'hello1' > "$3/in1"
    echo 'hello2' > "$3/in2"
    $shell "$2/test5_1.sh"
}

test5_2()
{
    local shell=$1
    echo 'hello1' > "$3/in1"
    echo 'hello2' > "$3/in2"
    $shell "$2/test5_2.sh"
}

test5_3()
{
    local shell=$1
    echo 'hello1' > "$3/in1"
    echo 'hello2' > "$3/in2"
    $shell "$2/test5_3.sh"
}

test6()
{
    local shell=$1
    $shell "$2/test6.sh"
}

test7_1()
{
    local shell=$1
    $shell "$2/test7_1.sh"
}

test7_2()
{
    local shell=$1
    $shell "$2/test7_2.sh"
}

test7_3()
{
    local shell=$1
    $shell "$2/test7_3.sh"
}

test8()
{
    local shell=$1
    $shell "$2/test8.sh"
}

test9_1()
{
    local shell=$1
    $shell "$2/test9_1.sh"
}

test9_2()
{
    local shell=$1
    $shell "$2/test9_2.sh"
}

test9_3()
{
    local shell=$1
    $shell "$2/test9_3.sh"
}

test_stdout()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > "$3/in1"
    $shell $2/test_stdout.sh
}

test_comments()
{
    local shell=$1
    $shell "$2/test_comments.sh"
}

test_function()
{
    local shell=$1
    $shell "$2/test_function.sh"
}

test_if()
{
    local shell=$1
    $shell $2/test_if.sh
}

test_if_2()
{
    local shell=$1
    $shell $2/test_if_2.sh
}

test_if_3()
{
    local shell=$1
    $shell $2/test_if_3.sh
}

test_cd()
{
    local shell=$1
    $shell $2/test_cd.sh
}

test_ifs()
{
    local shell=$1
    $shell $2/test_IFS.sh
}

test_loop()
{
    local shell=$1
    $shell $2/test_loop.sh
}

test_fd_1()
{
    local shell=$1
    $shell $2/test_fd_1.sh
}

test_commandline_args()
{
    local shell=$1
    $shell $2/test_commandline_args.sh a "b c" d
}

test_dynamic_exit()
{
    local shell=$1
    $shell $2/test_dynamic_exit.sh
}

test_break()
{
    local shell=$1
    $shell $2/test_break.sh
}

test_network_access_1()
{
    local shell=$1
    $shell $2/test_network_access_1.sh
}

test_network_access_2()
{
    local shell=$1
    $shell $2/test_network_access_2.sh
}

test_network_access_3()
{
    local shell=$1
    $shell $2/test_network_access_3.sh
}

test_local_vars_1()
{
    local shell=$1
    $shell $2/test_local_vars_1.sh
}

test_local_vars_2()
{
    local shell=$1
    $shell $2/test_local_vars_2.sh
}

test_local_vars_3()
{
    local shell=$1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in2
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in3
    $shell $2/test_local_vars_3.sh
}

test_local_vars_4()
{
    local shell=$1
    $shell $2/test_local_vars_4.sh
}

test_command_var_assignments_1(){
    local shell=$1
    $shell $2/test_command_var_assignments_1.sh
}

test_command_var_assignments_2(){
    local shell=$1
    $shell $2/test_command_var_assignments_2.sh
}

test_early_stop1()
{
    local shell=$1
    $shell $2/test_early_stop1.sh
}

test_early_stop2()
{
    local shell=$1
    $shell $2/test_early_stop2.sh
}

test_timed_execution()
{
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in2
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in3
    local shell=$1
    $shell $2/test_timed_execution.sh
}

test_timed_execution_2()
{
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in2
    local shell=$1
    $shell $2/test_timed_execution_2.sh
}

test_timed_loop()
{
    local shell=$1
    $shell $2/test_timed_loop.sh
}

test_scheduling_restart_1()
{
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in2
    local shell=$1
    $shell $2/test_scheduling_restart_1.sh
}

test_scheduling_restart_loop_1()
{
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in1
    echo $'foo\nbar\nbaz\nqux\nquux\nfoo\nbar' > $3/in_2_1
    local shell=$1
    $shell $2/test_scheduling_restart_loop_1.sh
}

## TODO: make more loop tests with nested loops and commands after the loop

# Arg parsing
if [ "$#" -gt 0 ]; then
    while [ "$#" -gt 0 ]; do
        case "$1" in
            test*)
                current_test=$1
                shift # Move past the test name
                execution_time_limit=""
                declare -a repetitions=()

                # Check if the next argument is numeric for the execution time limit
                if [[ "$1" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    execution_time_limit=$1
                    shift # Move past the execution time limit
                fi

                # All subsequent numeric arguments are repetitions
                while [[ "$1" =~ ^[0-9]+$ ]]; do
                    repetitions+=("$1")
                    shift # Move past each repetition
                done

                run_test "$current_test" "$execution_time_limit" "${repetitions[@]}"
                ;;
            *)
                echo "Unknown argument: $1"
                exit 1
                ;;
        esac
    done
else
    run_test test_single_command
    run_test test_local_vars_1
    run_test test_local_vars_2
    run_test test_local_vars_3
    run_test test_command_var_assignments_1
    run_test test_command_var_assignments_2
    # run_test test_scheduling_restart_1 5
    # run_test test_scheduling_restart_loop_1 5
    # run_test test_timed_execution 3
    # run_test test_timed_loop 5
    run_test test1_1 # "1 2 3 1" # 7
    run_test test1_2 #"1 2 2 1" # 6
    run_test test1_3 #"1 2 2 1" # 6
    run_test test2_1 #"1 1 1 1 1 1 1 1 1 1 1" # 10
    run_test test2_2 #"1 1 1 1 1 1 1 1 1 1 1" # 10
    run_test test2_3 #"1 1 1 1 1 1 1 1 1 1 1" # 10
    run_test test3_1 # "1 1 2 1 2" # 7
    run_test test3_2 # "1 1 2 1 2" # 7
    run_test test3_3 # "1 1 2 1 3" # 8
    run_test test4_1 #"1 2 1" # 4
    run_test test4_2 #"1 2 1" # 4
    run_test test4_3 #"1 2 1" # 4
    run_test test5_1 #"1 1 1" # 3
    run_test test5_2 #"1 1 1" # 3
    run_test test5_3 #"1 1 1" # 3
    # run_test test6
    run_test test7_1 #"1 1 1 1 1 1 1 1 1 1 1 1" # 12
    run_test test7_2 #"1 1 1 1 1 1 1 1 1 1 1 1" # 12
    run_test test7_3 #"1 1 1 1 1 1 1 1 1 1 1 1" # 12

    # for now we don't check for reps in tests 9_x
    run_test test9_1 # "1 2 1 1 1 2 2 2 2 2 1 1 1" # 19
    run_test test9_2 # "1 1 1 1 1 1 1 1 1 1 1 1 1" # 13
    run_test test9_3 # "1 1 1 1 1 1 1 1 2 2 1 1 1" # 15
    run_test test_stdout #"1 1 1 1 1 1" # 6
    run_test test_if
    run_test test_if_2
    run_test test_comments
    run_test test_function
    run_test test_loop
    run_test test_ifs
    run_test test_fd_1
    run_test test_commandline_args
    run_test test_break
    run_test test_network_access_1 #"1 2 2"
    run_test test_network_access_2 #"1 2 2 2"
    run_test test_network_access_3 #"1 2 2 2"
fi

if command -v lsb_release >/dev/null 2>&1; then
   distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ]; then
   distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# do different things depending on distro
case "$distro" in
    freebsd*)  
        # change sed to gsed
        sed () {
            gsed $@
        }
        ;;
    *)
        ;;
esac

echo -e "\n====================| Test Summary |====================\n"

echo "> Below follow the identical outputs:"
grep ": OK" "$output_dir/result_status" | awk '{print $1}' | sed 's/://g' | tee $output_dir/passed.log

echo "> Below follow the non-identical outputs:"
grep ": FAIL" "$output_dir/result_status" | awk '{print $1}' | sed 's/://g' | tee $output_dir/failed.log

TOTAL_TESTS=$(grep -cE "OK|FAIL" "$output_dir/result_status")
PASSED_TESTS=$(grep -c "OK" "$output_dir/result_status")

echo "========================================================"
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
echo "========================================================"
