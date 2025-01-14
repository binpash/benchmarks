import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def perform_pca_and_plot(dataframe):
    """
    Performs PCA on the numeric columns of the input dataframe and plots the first two principal components.

    Parameters:
        dataframe (pd.DataFrame): Input dataframe containing data for PCA.

    Returns:
        pd.DataFrame: A dataframe containing the principal components.
    """
    # Ensure numeric columns are selected for PCA
    numeric_cols = dataframe.select_dtypes(include=[np.number]).columns
    if numeric_cols.empty:
        raise ValueError("No numeric columns available in the dataframe for PCA.")
    print(f"Numeric columns selected for PCA: {numeric_cols}")

    # Drop rows with NaN values in numeric columns (if any)
    dataframe_numeric = dataframe[numeric_cols].dropna()

    # Standardize the data
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(dataframe_numeric)

    # Perform PCA
    pca = PCA(n_components=2)  # Reduce to 2 components for visualization
    principal_components = pca.fit_transform(data_scaled)

    # Create a new dataframe with the principal components
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    # Plot the results
    plt.figure(figsize=(10, 7))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.7)
    plt.title('PCA of Input DataFrame', fontsize=16)
    plt.xlabel('Principal Component 1', fontsize=12)
    plt.ylabel('Principal Component 2', fontsize=12)
    plt.grid(True)

    # Optionally, add labels for points (if 'benchmark' column exists)
    if 'benchmark' in dataframe.columns:
        for i, label in enumerate(dataframe['benchmark']):
            plt.annotate(label, (pca_df['PC1'][i], pca_df['PC2'][i]), fontsize=8, alpha=0.6)

    plt.savefig('pca_plot.pdf')

    return pca_df
