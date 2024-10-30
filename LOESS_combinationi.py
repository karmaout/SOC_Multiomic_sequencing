import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Load the LOESS regression data from both CSV files
pe_loess_file = "Region_PE_loess_regression_data.csv"  # Replace with the actual file path if needed
pr_loess_file = "Region_PR_loess_regression_data.csv"  # Replace with the actual file path if needed

pe_loess_data = pd.read_csv(pe_loess_file)
pr_loess_data = pd.read_csv(pr_loess_file)

# Create a new figure for the combined LOESS regression lines
fig, ax = plt.subplots(figsize=(12, 10))

# Plot LOESS regression line for PE data
ax.plot(pe_loess_data['x'], pe_loess_data['y'], color='blue', lw=2, linestyle='-', alpha=0.7, label='LOESS Regression (PE)', marker='o', markevery=50)

# Plot LOESS regression line for PR data
ax.plot(pr_loess_data['x'], pr_loess_data['y'], color='red', lw=2, linestyle='-', alpha=0.7, label='LOESS Regression (PR)', marker='s', markevery=50)

# Set the x-axis and y-axis range to match the previous plots (assuming the range was set from 0 to 350)
ax.set_xlim(0, 350)
ax.set_ylim(0, 350)

# Add x=y line for reference
ax.plot([0, 350], [0, 350], color='grey', linestyle='--', label='x = y')

# Customize plot labels and title with increased font sizes
ax.set_xlabel('Fos Positive ATAC', fontsize=16)
ax.set_ylabel('Fos Negative ATAC', fontsize=16)
ax.set_title('LOESS Regression Shift in Region', fontsize=20)

# Increase the font size of the tick labels
ax.tick_params(axis='both', which='major', labelsize=14)

# Add custom legend with increased font size
custom_lines = [
    Line2D([0], [0], color='blue', lw=2, linestyle='-', alpha=0.7, marker='o', label='LOESS Regression (PE)'),
    Line2D([0], [0], color='red', lw=2, linestyle='-', alpha=0.7, marker='s', label='LOESS Regression (PR)'),
    Line2D([0], [0], color='grey', lw=2, linestyle='--', label='x = y'),
]
ax.legend(handles=custom_lines, loc='upper left', fontsize=14)

# Save the plot as an SVG file for high-quality output
svg_file_path = "Region_Combined_LOESS_Regression_PE_PR.svg"  # Adjust the file name or path as needed
plt.savefig(svg_file_path, format='svg', dpi=300)  # Save as SVG with 300 DPI for better quality

# Show the plot
plt.tight_layout()
plt.show()
