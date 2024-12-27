#!/bin/bash

# Paths to files
REPORT_FILE=$(ls vps-audit-report-* 2>/dev/null | tail -n 1)
LOG_FILE=$(ls outputs/vps-audit-console.log 2>/dev/null | tail -n 1)
DIFF_FILE="log_report_diff.txt"

# Exit with a failure message
fail() {
    echo "Verification FAIL: $1"
    exit 1
}

# Ensure report and log files exist and are not empty
[[ -f "$REPORT_FILE" ]] || fail "No report file found. Please run the benchmark first."
[[ -f "$LOG_FILE" ]] || fail "No log file found in outputs directory. Check your setup."
[[ -s "$REPORT_FILE" ]] || fail "Report file is empty."
[[ -s "$LOG_FILE" ]] || fail "Log file is empty."

# Verify the presence of markers in the report
check_presence() {
    local type="$1"
    local file="$2"
    shift 2
    local items=("$@")

    for item in "${items[@]}"; do
        grep -qF "$item" "$file" || fail "Missing $type '$item' in $file."
    done
}

# Static markers
STATIC_MARKERS=(
    "VPS Security Audit Tool"
    "System Information"
    "Security Audit Results"
    "End of VPS Audit Report"
)

# Expected individual checks
EXPECTED_CHECKS=(
    "System Restart"
    "SSH Root Login"
    "SSH Password Auth"
    "SSH Port"
    "Firewall Status"
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

# Verify report structure and individual checks
check_presence "marker" "$REPORT_FILE" "${STATIC_MARKERS[@]}"
check_presence "check" "$REPORT_FILE" "${EXPECTED_CHECKS[@]}"

# Verify result markers ([PASS], [WARN], [FAIL])
RESULT_MARKERS=("[PASS]" "[WARN]" "[FAIL]")
RESULTS_COUNT=$(grep -Eo "$(IFS=\|; echo "${RESULT_MARKERS[*]}")" "$REPORT_FILE" | wc -l)

[[ "$RESULTS_COUNT" -gt 0 ]] || fail "No result markers ([PASS], [WARN], [FAIL]) found in the report."

# Disallowed markers
DISALLOWED_MARKERS=(
    "Error:"
    "command not found"
    "failed to retrieve"
    "No such file"
    "Cannot resolve"
)

# Check for disallowed markers
grep -F -f <(printf '%s\n' "${DISALLOWED_MARKERS[@]}") "$LOG_FILE" > "$DIFF_FILE" || true
[[ ! -s "$DIFF_FILE" ]] || fail "Found disallowed markers in log. Check '$DIFF_FILE' for details."

# All checks passed
echo "Verification PASS: All markers, checks, and comparisons are valid. No disallowed markers found."
exit 0
