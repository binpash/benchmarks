#!/bin/bash

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Get current timestamp for the report filename
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
REPORT_FILE="vps-audit-report.txt"

TEMP_DIR=$(mktemp -d)
export TEMP_DIR

# Cleanup on exit
trap "rm -rf $TEMP_DIR" EXIT

start_audit() {
    echo -e "${BLUE}${BOLD}VPS Security Audit Tool${NC}"
    echo -e "${GRAY}https://github.com/vernu/vps-audit${NC}"
    echo -e "${GRAY}Starting audit at $(date)${NC}\n"

    echo "VPS Security Audit Tool" > "$TEMP_DIR/001.txt"
    echo "https://github.com/vernu/vps-audit" >> "$TEMP_DIR/001.txt"
    echo "Starting audit at $(date)" >> "$TEMP_DIR/001.txt"
    echo "================================" >> "$TEMP_DIR/001.txt"
}

print_header_1() {
    echo -e "\n${BLUE}${BOLD}"System Information"${NC}"
    echo -e "\nSystem Information" >> "$TEMP_DIR/002.txt"
    echo "================================" >> "$TEMP_DIR/002.txt"
}


# print_info() {
#     local label="$1"
#     local value="$2"
#     echo -e "${BOLD}$label:${NC} $value"
#     echo "$label: $value"
# }


# Function to get Hostname
get_hostname() {
    echo "Hostname: $(hostname)" | tee "$TEMP_DIR/003.txt"
}

# Function to get OS information
get_os_info() {
    echo "Operating System: $(grep PRETTY_NAME /etc/os-release | cut -d'"' -f2)" | tee "$TEMP_DIR/004.txt"
}

# Function to get Kernel version
get_kernel_version() {
    echo "Kernel Version: $(uname -r)" | tee "$TEMP_DIR/005.txt"
}
# Function to get uptime
get_uptime() {
    echo "Uptime: $(uptime -p) (since $(uptime -s))" | tee "$TEMP_DIR/006.txt"
}

# Function to get CPU info
get_cpu_info() {
    echo "CPU Model: $(lscpu | grep 'Model name' | cut -d':' -f2 | xargs)" | tee "$TEMP_DIR/007.txt"
    echo "CPU Cores: $(nproc)" | tee -a "$TEMP_DIR/007.txt"
}


# Function to get memory info
get_memory_info() {
    echo "Total Memory: $(free -h | awk '/^Mem:/ {print $2}')" | tee "$TEMP_DIR/008.txt"
}

# Function to get disk space info
get_disk_info() {
    echo "Total Disk Space: $(df -h / | awk 'NR==2 {print $2}')" | tee "$TEMP_DIR/009.txt"
}

# Function to get public IP
get_public_ip() {
    echo "Public IP: $(curl -s https://api.ipify.org)" | tee "$TEMP_DIR/010.txt"
}

# Function to get load average
get_load_average() {
    echo "Load Average: $(uptime | awk -F'load average:' '{print $2}' | xargs)" | tee "$TEMP_DIR/011.txt"
    echo "" >> "$TEMP_DIR/011.txt"
}

# # Print system information
# print_info "Hostname" "$HOSTNAME"
# print_info "Operating System" "$OS_INFO"
# print_info "Kernel Version" "$KERNEL_VERSION"
# print_info "Uptime" "$UPTIME (since $UPTIME_SINCE)"
# print_info "CPU Model" "$CPU_INFO"
# print_info "CPU Cores" "$CPU_CORES"
# print_info "Total Memory" "$TOTAL_MEM"
# print_info "Total Disk Space" "$TOTAL_DISK"
# print_info "Public IP" "$PUBLIC_IP"
# print_info "Load Average" "$LOAD_AVERAGE"

# echo "" >> "$REPORT_FILE"

# Security Audit Section
print_header_2() {
    echo -e "\n${BLUE}${BOLD}Security Audit Results${NC}"
    echo -e "\nSecurity Audit Results" >> "$TEMP_DIR/012.txt"
    echo "================================" >> "$TEMP_DIR/012.txt"
}

check_security() {
    test_name="$1"
    status="$2"
    message="$3"
    outfile="$4"
    case $status in
        "PASS")
            echo -e "${GREEN}[PASS]${NC} $test_name ${GRAY}- $message${NC}\n" | tee "$TEMP_DIR/$4"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} $test_name ${GRAY}- $message${NC}\n" | tee "$TEMP_DIR/$4"
            ;;
        "FAIL")
            echo -e "${RED}[FAIL]${NC} $test_name ${GRAY}- $message${NC}\n" | tee "$TEMP_DIR/$4"
            ;;
    esac
}
export -f check_security
# Function to check system uptime
check_uptime() {
    UPTIME=$(uptime -p)
    UPTIME_SINCE=$(uptime -s)
    echo -e "\nSystem Uptime Information:" >> "$TEMP_DIR/013.txt"
    echo "Current uptime: $UPTIME" >> "$TEMP_DIR/013.txt"
    echo "System up since: $UPTIME_SINCE" >> "$TEMP_DIR/013.txt"
    echo "" >> "$TEMP_DIR/013.txt"
    echo -e "System Uptime: $UPTIME (since $UPTIME_SINCE)"
}

# Function to check if the system requires a restart
check_restart_required() {
    if [ -f /var/run/reboot-required ]; then
        check_security "System Restart" "WARN" "System requires a restart to apply updates" "014.txt"
    else
        check_security "System Restart" "PASS" "No restart required" "014.txt"
    fi
}

check_ssh_root_login() {
    SSH_CONFIG_OVERRIDES=$(grep "^Include" /etc/ssh/sshd_config 2>/dev/null | awk '{print $2}')

    if [ -n "$SSH_CONFIG_OVERRIDES" ] && [ -d "$(dirname "$SSH_CONFIG_OVERRIDES")" ]; then
        SSH_ROOT=$(grep "^PermitRootLogin" $SSH_CONFIG_OVERRIDES /etc/ssh/sshd_config 2>/dev/null | head -1 | awk '{print $2}')
    else
        SSH_ROOT=$(grep "^PermitRootLogin" /etc/ssh/sshd_config 2>/dev/null | head -1 | awk '{print $2}')
    fi
    if [ -z "$SSH_ROOT" ]; then
        SSH_ROOT="prohibit-password"
    fi
    if [ "$SSH_ROOT" = "no" ]; then
        check_security "SSH Root Login" "PASS" "Root login is properly disabled in SSH configuration" "015.txt"
    else
        check_security "SSH Root Login" "FAIL" "Root login is currently allowed - this is a security risk. Disable it in /etc/ssh/sshd_config" "015.txt"
    fi
}

# Function to check SSH password authentication
check_ssh_password_auth() {
    SSH_CONFIG_OVERRIDES=$(grep "^Include" /etc/ssh/sshd_config 2>/dev/null | awk '{print $2}')

    if [ -n "$SSH_CONFIG_OVERRIDES" ] && [ -d "$(dirname "$SSH_CONFIG_OVERRIDES")" ]; then
        SSH_PASSWORD=$(grep "^PasswordAuthentication" $SSH_CONFIG_OVERRIDES /etc/ssh/sshd_config 2>/dev/null | head -1 | awk '{print $2}')
    else
        SSH_PASSWORD=$(grep "^PasswordAuthentication" /etc/ssh/sshd_config 2>/dev/null | head -1 | awk '{print $2}')
    fi
    if [ -z "$SSH_PASSWORD" ]; then
        SSH_PASSWORD="yes"
    fi
    if [ "$SSH_PASSWORD" = "no" ]; then
        check_security "SSH Password Auth" "PASS" "Password authentication is disabled, key-based auth only" "016.txt"
    else
        check_security "SSH Password Auth" "FAIL" "Password authentication is enabled - consider using key-based authentication only" "016.txt"
    fi
}

# Function to check SSH port configuration
check_ssh_port() {
    UNPRIVILEGED_PORT_START=$(sysctl -n net.ipv4.ip_unprivileged_port_start)
    SSH_CONFIG_OVERRIDES=$(grep "^Include" /etc/ssh/sshd_config 2>/dev/null | awk '{print $2}')
    SSH_PORT=""
    if [ -n "$SSH_CONFIG_OVERRIDES" ] && [ -d "$(dirname "$SSH_CONFIG_OVERRIDES")" ]; then
        SSH_PORT=$(grep "^Port" $SSH_CONFIG_OVERRIDES /etc/ssh/sshd_config 2>/dev/null | head -1 | awk '{print $2}')
    else
        SSH_PORT=$(grep "^Port" /etc/ssh/sshd_config 2>/dev/null | head -1 | awk '{print $2}')
    fi
    if [ -z "$SSH_PORT" ]; then
        SSH_PORT="22"
    fi
    if [ "$SSH_PORT" = "22" ]; then
        check_security "SSH Port" "WARN" "Using default port 22 - consider changing to a non-standard port for security by obscurity" "017.txt"
    elif [ "$SSH_PORT" -ge "$UNPRIVILEGED_PORT_START" ]; then
        check_security "SSH Port" "FAIL" "Using unprivileged port $SSH_PORT - use a port below $UNPRIVILEGED_PORT_START for better security" "017.txt"
    else
        check_security "SSH Port" "PASS" "Using non-default port $SSH_PORT which helps prevent automated attacks" "017.txt"
    fi
}


# Function to check unattended upgrades
check_unattended_upgrades() {
    if dpkg -l | grep -q "unattended-upgrades"; then
        check_security "Unattended Upgrades" "PASS" "Automatic security updates are configured" "019.txt"
    else
        check_security "Unattended Upgrades" "FAIL" "Automatic security updates are not configured - system may miss critical updates" "019.txt"
    fi
}

# Function to check fail2ban
check_fail2ban() {
    if dpkg -l | grep -q "fail2ban"; then
        if pgrep -x "fail2ban-server" >/dev/null 2>&1; then
            check_security "Fail2ban" "PASS" "Brute force protection is active and running" "020.txt"
        else
            check_security "Fail2ban" "WARN" "Fail2ban is installed but not running - brute force protection is disabled" "020.txt"
        fi
    else
        check_security "Fail2ban" "FAIL" "No brute force protection installed - system is vulnerable to login attacks" "020.txt"
    fi
}

# Function to check failed login attempts
check_failed_logins() {
    failed_logins=$(grep "Failed password" /var/log/auth.log 2>/dev/null | wc -l)
    if [ "$failed_logins" -lt 10 ]; then
        check_security "Failed Logins" "PASS" "Only $failed_logins failed login attempts detected - this is within normal range" "021.txt"
    elif [ "$failed_logins" -lt 50 ]; then
        check_security "Failed Logins" "WARN" "$failed_logins failed login attempts detected - might indicate breach attempts" "021.txt"
    else
        check_security "Failed Logins" "FAIL" "$failed_logins failed login attempts detected - possible brute force attack in progress" "021.txt"
    fi
}

# Function to check system updates
check_system_updates() {
    updates=$(apt-get -s upgrade 2>/dev/null | grep -P '^\d+ upgraded' | cut -d" " -f1)
    if [ "$updates" -eq 0 ]; then
        check_security "System Updates" "PASS" "All system packages are up to date" "022.txt"
    else
        check_security "System Updates" "FAIL" "$updates security updates available - system is vulnerable to known exploits" "022.txt"
    fi
}
# Function to check running services
check_running_services() {
    services=$(ps --no-headers -eo cmd | wc -l)

    if [ "$services" -lt 20 ]; then
        check_security "Running Services" "PASS" "Running minimal services ($services) - good for security" "023.txt"
    elif [ "$services" -lt 40 ]; then
        check_security "Running Services" "WARN" "$services services running - consider reducing attack surface" "023.txt"
    else
        check_security "Running Services" "FAIL" "Too many services running ($services) - increases attack surface" "023.txt"
    fi
}

# Function to check open ports
check_ports() {
    if command -v netstat >/dev/null 2>&1; then
        listening_ports=$(netstat -tuln | grep LISTEN | awk '{print $4}')
    elif command -v ss >/dev/null 2>&1; then
        listening_ports=$(ss -tuln | grep LISTEN | awk '{print $5}')
    else
        check_security "Port Scanning" "FAIL" "Neither 'netstat' nor 'ss' is available on this system." "024.txt"
        return
    fi

    if [ -n "$listening_ports" ]; then
        public_ports=$(echo "$listening_ports" | awk -F':' '{print $NF}' | sort -n | uniq | tr '\n' ',' | sed 's/,$//')
        port_count=$(echo "$public_ports" | tr ',' '\n' | wc -w)
        internet_ports=$(echo "$public_ports" | tr ',' '\n' | wc -w)

        if [ "$port_count" -lt 10 ] && [ "$internet_ports" -lt 3 ]; then
            check_security "Port Security" "PASS" "Good configuration (Total: $port_count, Public: $internet_ports accessible ports): $public_ports" "025.txt"
        elif [ "$port_count" -lt 20 ] && [ "$internet_ports" -lt 5 ]; then
            check_security "Port Security" "WARN" "Review recommended (Total: $port_count, Public: $internet_ports accessible ports): $public_ports" "025.txt"
        else
            check_security "Port Security" "FAIL" "High exposure (Total: $port_count, Public: $internet_ports accessible ports): $public_ports" "025.txt"
        fi
    else
        check_security "Port Scanning" "WARN" "Port scanning failed due to missing tools. Ensure 'ss' or 'netstat' is installed." "025.txt"
    fi
}

# Function to check disk usage
check_disk_usage() {
    disk_total=$(df -h / | awk 'NR==2 {print $2}')
    disk_used=$(df -h / | awk 'NR==2 {print $3}')
    disk_avail=$(df -h / | awk 'NR==2 {print $4}')
    disk_usage=$(df -h / | awk 'NR==2 {print int($5)}')

    if [ "$disk_usage" -lt 50 ]; then
        check_security "Disk Usage" "PASS" "Healthy disk space available (${disk_usage}% used - Used: ${disk_used} of ${disk_total}, Available: ${disk_avail})" "026.txt"
    elif [ "$disk_usage" -lt 80 ]; then
        check_security "Disk Usage" "WARN" "Disk space usage is moderate (${disk_usage}% used - Used: ${disk_used} of ${disk_total}, Available: ${disk_avail})" "026.txt"
    else
        check_security "Disk Usage" "FAIL" "Critical disk space usage (${disk_usage}% used - Used: ${disk_used} of ${disk_total}, Available: ${disk_avail})" "026.txt"
    fi
}

check_memory_usage() {
    MEM_TOTAL=$(free -h | awk '/^Mem:/ {print $2}')
    MEM_USED=$(free -h | awk '/^Mem:/ {print $3}')
    MEM_AVAIL=$(free -h | awk '/^Mem:/ {print $7}')
    MEM_USAGE=$(free | awk '/^Mem:/ {printf "%.0f", $3/$2 * 100}')
    if [ "$MEM_USAGE" -lt 50 ]; then
        check_security "Memory Usage" "PASS" "Healthy memory usage (${MEM_USAGE}% used - Used: ${MEM_USED} of ${MEM_TOTAL}, Available: ${MEM_AVAIL})" "027.txt"
    elif [ "$MEM_USAGE" -lt 80 ]; then
        check_security "Memory Usage" "WARN" "Moderate memory usage (${MEM_USAGE}% used - Used: ${MEM_USED} of ${MEM_TOTAL}, Available: ${MEM_AVAIL})" "027.txt"
    else
        check_security "Memory Usage" "FAIL" "Critical memory usage (${MEM_USAGE}% used - Used: ${MEM_USED} of ${MEM_TOTAL}, Available: ${MEM_AVAIL})" "027.txt"
    fi
}


# Function to check CPU usage
check_cpu_usage() {
    cpu_cores=$(nproc)
    cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print int($2)}')
    cpu_idle=$(top -bn1 | grep "Cpu(s)" | awk '{print int($8)}')
    cpu_load=$(uptime | awk -F'load average:' '{ print $2 }' | awk -F',' '{ print $1 }' | tr -d ' ')
    
    if [ "$cpu_usage" -lt 50 ]; then
        check_security "CPU Usage" "PASS" "Healthy CPU usage (${cpu_usage}% used - Active: ${cpu_usage}%, Idle: ${cpu_idle}%, Load: ${cpu_load}, Cores: ${cpu_cores})" "028.txt"
    elif [ "$cpu_usage" -lt 80 ]; then
        check_security "CPU Usage" "WARN" "Moderate CPU usage (${cpu_usage}% used - Active: ${cpu_usage}%, Idle: ${cpu_idle}%, Load: ${cpu_load}, Cores: ${cpu_cores})" "028.txt"
    else
        check_security "CPU Usage" "FAIL" "Critical CPU usage (${cpu_usage}% used - Active: ${cpu_usage}%, Idle: ${cpu_idle}%, Load: ${cpu_load}, Cores: ${cpu_cores})" "028.txt"
    fi
}

# Function to check password policy
check_password_policy() {
    if [ -f "/etc/security/pwquality.conf" ]; then
        if grep -q "minlen.*12" /etc/security/pwquality.conf; then
            check_security "Password Policy" "PASS" "Strong password policy is enforced" "030.txt"
        else
            check_security "Password Policy" "FAIL" "Weak password policy - passwords may be too simple" "030.txt"
        fi
    else
        check_security "Password Policy" "FAIL" "No password policy configured - system accepts weak passwords" "030.txt"
    fi
}

# Function to check for suspicious SUID files
check_suid_files() {
    local suid_files
    COMMON_SUID_PATHS='^/usr/bin/|^/bin/|^/sbin/|^/usr/sbin/|^/usr/lib|^/usr/libexec'
    KNOWN_SUID_BINS='ping$|sudo$|mount$|umount$|su$|passwd$|chsh$|newgrp$|gpasswd$|chfn$'

    suid_files=$(find / -type f -perm -4000 2>/dev/null | \
        grep -v -E "$COMMON_SUID_PATHS" | \
        grep -v -E "$KNOWN_SUID_BINS" | \
        wc -l)

    if [ "$suid_files" -eq 0 ]; then
        check_security "SUID Files" "PASS" "No suspicious SUID files found - good security practice" "031.txt"
    else
        check_security "SUID Files" "WARN" "Found $suid_files SUID files outside standard locations - verify if legitimate" "031.txt"
    fi
}

# Call the final report generation function
# generate_final_report

# # Notify user of completion
# echo -e "\n${GREEN}${BOLD}VPS audit complete.${NC} Full report saved to: $REPORT_FILE"
# echo -e "Review $REPORT_FILE for detailed recommendations."

# Add system information summary to report
generate_sysinfo() {
    echo "================================" >> "$TEMP_DIR/032.txt"
    echo "System Information Summary:" >> "$TEMP_DIR/032.txt"
    echo "Hostname: $(hostname)" >> "$TEMP_DIR/032.txt"
    echo "Kernel: $(uname -r)" >> "$TEMP_DIR/032.txt"
    echo "OS: $(cat /etc/os-release | grep PRETTY_NAME | cut -d'"' -f2)" >> "$TEMP_DIR/032.txt"
    echo "CPU Cores: $(nproc)" >> "$TEMP_DIR/032.txt"
    echo "Total Memory: $(free -h | awk '/^Mem:/ {print $2}')" >> "$TEMP_DIR/032.txt"
    echo "Total Disk Space: $(df -h / | awk 'NR==2 {print $2}')" >> "$TEMP_DIR/032.txt"
    echo "================================" >> "$TEMP_DIR/032.txt"
    echo "================================" >> "$TEMP_DIR/032.txt"
    echo "End of VPS Audit Report" >> "$TEMP_DIR/032.txt"
    echo "Please review all failed checks and implement the recommended fixes." >> "$TEMP_DIR/032.txt"
}

export -f start_audit get_os_info get_kernel_version get_hostname get_uptime \
          get_cpu_info get_memory_info get_disk_info get_public_ip \
          get_load_average check_uptime check_restart_required \
          check_ssh_root_login check_ssh_password_auth check_ssh_port \
          check_firewall_status check_unattended_upgrades check_fail2ban \
          check_failed_logins check_system_updates \
          check_running_services check_ports check_disk_usage \
          check_memory_usage check_cpu_usage check_sudo_logging \
          check_password_policy check_suid_files \
          generate_sysinfo print_header_1 print_header_2

parallel ::: print_header_1 start_audit  \
             get_os_info get_kernel_version get_hostname get_uptime \
             get_cpu_info get_memory_info get_disk_info get_public_ip \
             get_load_average print_header_2 check_uptime check_restart_required \
             check_ssh_root_login check_ssh_password_auth check_ssh_port \
             check_firewall_status check_unattended_upgrades check_fail2ban \
             check_failed_logins check_system_updates \
             check_running_services check_ports check_disk_usage \
             check_memory_usage check_cpu_usage check_sudo_logging \
             check_password_policy check_suid_files generate_sysinfo

cat "$TEMP_DIR"/*.txt >> "$REPORT_FILE"

echo -e "\n${GREEN}${BOLD}System Information collected.${NC}"
echo -e "\nVPS audit complete. Full report saved to $REPORT_FILE"
echo -e "Review $REPORT_FILE for detailed recommendations."
