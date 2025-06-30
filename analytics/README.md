## analytics

This benchmark analyzes web server logs and network capture files to extract information about client behavior, error patterns, and protocol-level content.

### Inputs

- `inputs/nginx-logs/`: A directory containing NGINX access log files.
- `inputs/pcaps/`: A directory containing pcap (packet capture) files.

### Running

Two analysis scripts are executed:

1. **NGINX Log Analysis**:
   - Each access log is processed to extract metrics such as:
     - Frequency of HTTP response codes (e.g., 404, 502)
     - URLs requested
     - IPs requesting specific paths
   - The output is saved to `outputs/nginx_<size>/<log_filename>`.

2. **PCAP Analysis**:
   - Each pcap file is inspected using `tcpdump` to extract:
     - DNS queries
     - HTTP requests and hosts
     - Password patterns in plaintext traffic
     - Telnet login attempts
   - The output is saved to `outputs/pcaps_<size>/<filename>.log`.

3. **IP-ASN Annotation and ASN Popularity Analysis:**:
   - Each log file containing JSON-structured IP information is:
     - Annotated using zannotate with BGP routing information from a provided MRT file.
     - Parsed using jq to extract:
      - Source IP addresses
      - Associated Autonomous System Numbers (ASNs)
   - These are joined and processed to compute the popularity of ASNs, based on how frequently they appear in the annotated logs.
   - The output is saved to 
    ```bash
   outputs/asn_<size>/<log_filename>/
    ├── annotated.json         # Full annotation output
    ├── ips.txt                # Extracted IPs
    ├── asns.txt               # Corresponding ASNs
    └── as_popularity.csv      # ASN popularity (ranked by frequency)
    ```
    
### Validation

Correctness is determined by computing the MD5 hash of each output file and comparing it against a reference hash stored in `hashes/`.

### References

- https://zenodo.org/records/4743746
- https://www.usenix.org/system/files/atc20-raghavan.pdf
