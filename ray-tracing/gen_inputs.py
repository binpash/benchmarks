import os
import random
import csv

output_dir = "ray-tracing/inputs"
os.makedirs(output_dir, exist_ok=True)

headers = ["pathID", "timestamp", "hop", "value"]

info1_path = os.path.join(output_dir, "1.INFO")
with open(info1_path, "w") as f:
    f.write("[RAY] pathID,timestamp,hop,value\n")
    for i in range(10000):
        f.write(f"[RAY] {20613314 + i},{random.randint(1000,5000)},{random.randint(1,20)},{random.random():.4f}\n")

for idx in range(2, 6):
    path = os.path.join(output_dir, f"{idx}.INFO")
    with open(path, "w") as f:
        for _ in range(10000):
            f.write(f"[RAY] {20613314 + random.randint(0,4)},{random.randint(1000,5000)},{random.randint(1,20)},{random.random():.4f}\n")

output_dir
