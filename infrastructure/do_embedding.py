import json
import glob
import os
from openai import OpenAI
import pandas as pd
import dotenv

# Load environment variables
dotenv.load_dotenv()

client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# Set your OpenAI API key

def read_json_and_generate_embeddings(json_file):
    """
    Reads the JSON file, combines script contents for each benchmark, sends them to the OpenAI embedding API,
    and returns a dictionary of dataframes containing embeddings for each benchmark.
    """
    # Read the JSON file
    with open(json_file, 'r') as file:
        data = json.load(file)

    embedding_df = pd.DataFrame(columns=["benchmark", "embedding"])

    # Process each benchmark
    for benchmark, details in data.items():
        print(f"Processing benchmark: {benchmark}")

        # Combine all script contents into a single string
        scripts_globs = details.get("scripts", [])
        combined_script = ""
        for script_glob in scripts_globs:
            for script_file in glob.glob(f"../{script_glob}"):
                # Skip auxiliary files
                if "test" in script_file:
                    continue
                if "negate" in script_file:
                    continue
                if "header" in script_file:
                    continue
                with open(script_file, 'r') as f:
                    combined_script += f.read() + "\n"  # Append content

        print(len(combined_script))
        # Generate embedding using OpenAI's API
        try:
            response = client.embeddings.create(model="text-embedding-ada-002",  # Use a suitable model for embedding
            input=combined_script)
            embedding = response.data[0].embedding
        except Exception as e:
            print(f"Error generating embedding for {benchmark}: {e}")
            continue

        # Create a dataframe to hold the benchmark and its embedding
        embedding_df = embedding_df._append({"benchmark": benchmark, "embedding": embedding}, ignore_index=True)

    return embedding_df

# Example usage
if __name__ == "__main__":
    json_file = "./data/script-globs.json"
    embeddings_df = read_json_and_generate_embeddings(json_file)

    # Save or inspect the results
    print(embeddings_df)
    embeddings_df.to_csv("./data/embeddings.csv", index=False)
