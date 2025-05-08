## vps-audit

This benchmark performs an automated security assessment of a Linux-based VPS
(Virtual Private Server), focusing on system configurations, services, resource
usage, and known best practices.

### Overview

The audit runs two scripts:
- `vps-audit.sh`: Performs the main security and configuration audit.
- `vps-audit-negate.sh`: Performs a similar audit with flipped control logic or edge-case validation.

### Features Checked

- System uptime and restart requirements
- SSH configuration (root login, password authentication, port)
- Firewall configuration (commented out by default)
- Automatic updates (`unattended-upgrades`)
- Brute-force protection (`fail2ban`)
- Failed login attempts
- Installed security updates
- Running services and open ports
- Disk, memory, and CPU usage
- Password policy and suspicious SUID files

### Output

Each script produces a report file:
- `vps-audit-report.txt`
- `vps-audit-negate-report.txt`

### Validation

The `validate.py`:
- Canonicalizes output (removing timestamps, formatting noise)
- Compares SHA-256 hashes against reference files in `hash/`

### References

- https://github.com/vernu/vps-audit
