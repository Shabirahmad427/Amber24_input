import pandas as pd
import matplotlib.pyplot as plt

# === CONFIGURATION ===
input_file = 'output_hbonding.dat'
top_n = 10
save_fig = True
output_fig = 'top_hbonds.png'

# === DATA PARSING ===
donors = []
acceptors = []
occupancies = []

with open(input_file, 'r') as file:
    for line in file:
        if '|' in line and not line.startswith('Donor'):
            parts = line.strip().split('|')
            left = parts[0].strip()
            right = parts[1].strip()
            if len(left.split()) >= 2:
                donor, acceptor = left.split(None, 1)
                donors.append(donor.strip())
                acceptors.append(acceptor.strip())
                occupancies.append(float(right))

# Create DataFrame
df = pd.DataFrame({
    'Pair': [f'{d} → {a}' for d, a in zip(donors, acceptors)],
    'Occupancy': occupancies
})

# Sort and select top N
df_sorted = df.sort_values(by='Occupancy', ascending=False).head(top_n)

# === PLOT SETUP ===
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'Arial',  # Set Arial font
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11
})

fig, ax = plt.subplots(figsize=(10, 6))

bars = ax.barh(df_sorted['Pair'], df_sorted['Occupancy'],
               color='skyblue', edgecolor='black')

ax.set_xlabel('Occupancy (%)')
ax.set_title(f'')
ax.invert_yaxis()
ax.grid(axis='x', linestyle='--', alpha=0.6)

# Add text labels to bars
for bar in bars:
    width = bar.get_width()
    ax.text(width + 0.2, bar.get_y() + bar.get_height()/2,
            f'{width:.2f}', va='center', fontsize=10)

plt.tight_layout()

# Save or show
if save_fig:
    plt.savefig(output_fig, dpi=300, bbox_inches='tight')
plt.show()
