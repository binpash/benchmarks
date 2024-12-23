#!/bin/bash

# Paths to files
REPORT_FILE=$(ls vps-audit-report-* 2>/dev/null | tail -n 1)
LOG_FILE=$(ls outputs/vps-audit-console.log 2>/dev/null | tail -n 1)
DIFF_FILE="log_report_diff.txt"

# Ensure report and log files exist
if [[ ! -f "$REPORT_FILE" ]]; then
    echo "Verification FAIL: No report file found. Please run the benchmark first."
    exit 1
fi

if [[ ! -f "$LOG_FILE" ]]; then
    echo "Verification FAIL: No log file found in outputs directory. Check your setup."
    exit 1
fi

# Ensure the report file is not empty
if [[ ! -s "$REPORT_FILE" ]]; then
    echo "Verification FAIL: Report file is empty."
    exit 1
fi

# Ensure the log file is not empty
if [[ ! -s "$LOG_FILE" ]]; then
    echo "Verification FAIL: Log file is empty."
    exit 1
fi

# Define static markers to verify the report structure
STATIC_MARKERS=(
    "VPS Security Audit Tool"
    "System Information"
    "Security Audit Results"
    "End of VPS Audit Report"
)

# Define individual checks that must appear in the report
EXPECTED_CHECKS=(
    "System Restart"
    "SSH Root Login"
    "SSH Password Auth"
    "SSH Port"
    "Firewall Status"
    "Intrusion Prevention"
    "Failed Logins"
    "System Updates"
    "Port Security"
    "Disk Usage"
    "Memory Usage"
    "CPU Usage"
    "Sudo Logging"
    "Password Policy"
    "SUID Files"
)

# Define disallowed markers
DISALLOWED_MARKERS=(
    "Error:"
    "command not found"
    "failed to retrieve"
    "No such file" 
    "Cannot resolve"
    # "Permission denied" # may have to be removed, else we need sudo
)

# Verify the presence of static markers
for marker in "${STATIC_MARKERS[@]}"; do
    if ! grep -qF "$marker" "$REPORT_FILE"; then
        echo "Verification FAIL: Missing marker '$marker' in report."
        exit 1
    fi
done

# Verify the presence of individual checks
for check in "${EXPECTED_CHECKS[@]}"; do
    if ! grep -qF "$check" "$REPORT_FILE"; then
        echo "Verification FAIL: Missing check '$check' in report."
        exit 1
    fi
done

# Verify the presence of result markers ([PASS], [WARN], [FAIL])
EXPECTED_RESULTS=("[PASS]" "[WARN]" "[FAIL]")
RESULTS_COUNT=0

for result in "${EXPECTED_RESULTS[@]}"; do
    COUNT=$(grep -o "$result" "$REPORT_FILE" | wc -l)
    if [[ "$COUNT" -gt 0 ]]; then
        RESULTS_COUNT=$((RESULTS_COUNT + COUNT))
    fi
done

if [[ "$RESULTS_COUNT" -eq 0 ]]; then
    echo "Verification FAIL: No result markers ([PASS], [WARN], [FAIL]) found in the report."
    exit 1
fi

# Check for disallowed markers in the diff
grep -F -f <(printf '%s\n' "${DISALLOWED_MARKERS[@]}") "$LOG_FILE" > "$DIFF_FILE"

if [[ -s "$DIFF_FILE" ]]; then
    echo "Verification FAIL: Found disallowed markers in log. Check '$DIFF_FILE' for details."
    exit 1
fi

# If all checks pass
echo "Verification PASS: All markers, checks, and comparisons are valid. No disallowed markers found."
exit 0
