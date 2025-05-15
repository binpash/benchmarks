import os
import random
import argparse

def generate_info_files(output_dir: str = "inputs", num_lines: int = 10000, seed: int = 42) -> str:
    os.makedirs(output_dir, exist_ok=True)
    random.seed(42)

    info1_path = os.path.join(output_dir, "1.INFO")
    with open(info1_path, "w") as f:
        f.write("[RAY] pathID,timestamp,hop,value\n")
        for i in range(num_lines):
            f.write(f"[RAY] {20613314 + i},{random.randint(1000, 5000)},{random.randint(1, 20)},{random.random():.4f}\n")

    for idx in range(2, 6):
        path = os.path.join(output_dir, f"{idx}.INFO")
        with open(path, "w") as f:
            for _ in range(num_lines):
                f.write(f"[RAY] {20613314 + random.randint(0, 4)},{random.randint(1000, 5000)},{random.randint(1, 20)},{random.random():.4f}\n")

    return output_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate deterministic .INFO files for RAY analysis.")
    parser.add_argument("output_dir", nargs="?", default="inputs", help="Directory to store generated .INFO files")
    parser.add_argument("num_lines", nargs="?", type=int, default=10000, help="Number of lines per .INFO file")

    args = parser.parse_args()

    result_dir = generate_info_files(args.output_dir, args.num_lines)
    print(f"Output directory: {result_dir}")
