#!/usr/bin/env python3

import pandas as pd
from do_pca import perform_pca_and_plot
import ast
from project_root import get_project_root
import sys
from pathlib import Path

OUTPUT = sys.argv[1] if len(sys.argv) > 1 else None

# make OUTPUT absolute path
if OUTPUT:
    OUTPUT = Path(OUTPUT).absolute()
else:
    OUTPUT = Path(__file__).parent / 'pca-row-plot.pdf'

root = get_project_root()


embedding_df = pd.read_csv(root / 'infrastructure/data/embeddings.csv')
embedding_df['embedding'] = embedding_df['embedding'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

# Embedding is a list of numbers, turn them into columns
embedding_df = pd.concat([embedding_df['benchmark'], embedding_df['embedding'].apply(pd.Series)], axis=1)

big_bench = embedding_df

perform_pca_and_plot(big_bench, embedding_df, str(OUTPUT))
