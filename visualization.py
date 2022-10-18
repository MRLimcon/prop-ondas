import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("./results/data.csv", index_col=0)
time = []
distances = {}
values = []

fig, axes = plt.subplots(nrows=6, ncols=1)
for column in data.columns:
    if column == "Time (s)":
        time = data[column]
    else:
        distances[column] = data[column]
        values.append(column)

for i, distance in enumerate(values):
    axes[i].plot(time, distances[distance])
    axes[i].set_title(distance)

fig.tight_layout()
plt.show()