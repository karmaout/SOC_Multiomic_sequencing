import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from matplotlib.lines import Line2D
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess  # Import LOWESS from statsmodels
from adjustText import adjust_text  # Import the adjustText library

# Load the data from the file
file_path = "REGION_PE_anno-sorted4_foldchange.txt"  # Replace with the actual path to your file
data = pd.read_csv(file_path, sep='\t')

# Ensure relevant columns are converted to the correct data types
data['fos_pos_atac'] = data['fos_pos_atac'].astype(float)
data['fos_neg_atac'] = data['fos_neg_atac'].astype(float)
data['motif_weight'] = data['motif_weight'].astype(float)

# Calculate deviation as the absolute difference between fos_pos_atac and fos_neg_atac
data['deviation'] = abs(data['fos_pos_atac'] - data['fos_neg_atac'])

# Extract relevant columns for plotting
fos_pos_atac = data['fos_pos_atac']
fos_neg_atac = data['fos_neg_atac']
motif_weight = data['motif_weight']

# Create a colormap based on motif_weight, from yellow to red
norm = mcolors.Normalize(vmin=motif_weight.min(), vmax=motif_weight.max())
cmap = plt.cm.YlOrRd

# Create the figure and the scatter plot
fig, ax = plt.subplots(figsize=(12, 10))  # Increase the figure size for better readability

# Plot all points in light grey
ax.scatter(fos_pos_atac, fos_neg_atac, color='lightgrey', marker='o', s=20, label='All ATAC Peaks', alpha=0.5)

# Highlight points based on motif_weight in yellow to red gradient
scatter = ax.scatter(fos_pos_atac, fos_neg_atac, c=motif_weight, cmap=cmap, s=motif_weight, alpha=0.8, label='Motif Weight')

# Set the x-axis and y-axis range
ax.set_xlim(0, 400)
ax.set_ylim(0, 400)

# Extend the x=y line across the full plot
ax.plot([0, 400], [0, 400], 'grey', linestyle='--', label='x = y')

# Fit a LOWESS regression model using statsmodels
lowess_result = lowess(endog=fos_neg_atac, exog=fos_pos_atac, frac=0.3)
sorted_indices = np.argsort(lowess_result[:, 0])
sorted_lowess = lowess_result[sorted_indices]

# Plot the original LOESS regression line
ax.plot(sorted_lowess[:, 0], sorted_lowess[:, 1], color='blue', lw=2, linestyle='-', label='LOESS Regression (Observed)')

# Extract slope and points at boundaries
x_start, y_start = sorted_lowess[0]
x_end, y_end = sorted_lowess[-1]
slope_start = (sorted_lowess[1, 1] - sorted_lowess[0, 1]) / (sorted_lowess[1, 0] - sorted_lowess[0, 0])
slope_end = (sorted_lowess[-1, 1] - sorted_lowess[-2, 1]) / (sorted_lowess[-1, 0] - sorted_lowess[-2, 0])

# Create extended points based on slope
extended_x_start = np.linspace(0, x_start, 100)
extended_y_start = y_start + slope_start * (extended_x_start - x_start)

extended_x_end = np.linspace(x_end, 400, 100)
extended_y_end = y_end + slope_end * (extended_x_end - x_end)

# Plot the extended regression lines with dashed lines
ax.plot(extended_x_start, extended_y_start, color='blue', lw=2, linestyle='--', label='LOESS Regression (Extended)')
ax.plot(extended_x_end, extended_y_end, color='blue', lw=2, linestyle='--')

# Save the extended LOWESS regression line data as a CSV file for further use
extended_loess_data = pd.DataFrame({'x': np.concatenate([extended_x_start, sorted_lowess[:, 0], extended_x_end]),
                                    'y': np.concatenate([extended_y_start, sorted_lowess[:, 1], extended_y_end])})
extended_loess_data.to_csv('REGION_PE_loess_regression_data.csv', index=False)

# Add color bar for motif_weight
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Fos Motif Weight', fontsize=22)

# Customize plot labels and title
ax.set_xlabel('Fos Positive ATAC', fontsize=22)
ax.set_ylabel('Fos Negative ATAC', fontsize=22)
ax.set_title('REGION PE Scatter Plot of Fos Positive and Negative ATAC Peaks', fontsize=18)

# Add a custom legend
custom_lines = [
    Line2D([0], [0], color='blue', lw=2, linestyle='-', label='LOESS Regression (Observed)'),
    Line2D([0], [0], color='blue', lw=2, linestyle='--', label='LOESS Regression (Extended)'),
    Line2D([0], [0], color='grey', lw=2, linestyle='--', label='x = y'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgrey', markersize=10, label='All ATAC Peaks'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Motif Weight')
]
ax.legend(handles=custom_lines, loc='upper left', fontsize=14, frameon=True)

# Save the plot as a vector file (SVG)
svg_file_path = 'REGION_PE_fos_atac_plot_no_label.tiff'  # Adjust the path to where you want to save the file
plt.savefig(svg_file_path, format='tiff')  # Save as TIFF format

# Show the plot
plt.tight_layout()
plt.show()
