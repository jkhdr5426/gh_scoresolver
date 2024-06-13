import matplotlib.pyplot as plt

file_path = '/Users/edwardkim/Desktop/gh_scoresolver/pdb_6/hist/pmf.txt'
coordinates = []
free_energies = []
probabilities = []
with open(file_path, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        columns = line.split()
        if len(columns) == 5:
            coor = float(columns[0])
            free = float(columns[1]) if columns[1] != 'nan' else None
            prob = float(columns[3]) if columns[3] != 'nan' else None
            coordinates.append(coor)
            free_energies.append(free)
            probabilities.append(prob)
fig, ax1 = plt.subplots()
ax1.set_xlabel('Coordinate')
ax1.set_ylabel('Free Energy', color='tab:blue')
ax1.plot(coordinates, free_energies, 'b-', label='Free Energy')
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax2 = ax1.twinx()
ax2.set_ylabel('Probability', color='tab:red')
ax2.plot(coordinates, probabilities, 'r-', label='Probability')
ax2.tick_params(axis='y', labelcolor='tab:red')
fig.tight_layout()
fig.legend(loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
plt.title('Histogram Data Plot')
plt.show()
