import re
import json
import sys
import random
seed = 42
random.seed(seed)
# Updated regex: allow "-" or number for response size
pattern = re.compile(
    r'(?P<ip>\S+) '                       # IP
    r'\S+ \S+ '                           # identd and user
    r'\[(?P<time>[^\]]+)\] '              # time
    r'"(?P<method>\S+) (?P<path>.*?) (?P<protocol>[^"]+)" '  # request
    r'(?P<status>\d{3}) '                 # status
    r'(?P<size>\d+|-)'                    # size can be a number or "-"
    r' "(?P<referer>[^"]*)" '             # referer
    r'"(?P<user_agent>[^"]*)"'            # user agent
)

def parse_line(line):
    match = pattern.match(line)
    if match:
        data = match.groupdict()

        # Add default values for optional fields
        for field in ["referer", "user_agent"]:
            if not data.get(field) or data[field] in ["-", "null", ""]:
                data[field] = "UNKNOWN"

        # Add a randomly generated ASN (for example in the public 32-bit range)
        random_asn = random.randint(131072, 4199999999)
        data["asn"] = str(random_asn)

        return data
    else:
        return None


def main(file_path):
    with open(file_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            data = parse_line(line)
            if data:
                print(json.dumps(data, ensure_ascii=False))
            else:
                continue


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python parse_nginx.py <access.log>")
    else:
        main(sys.argv[1])
