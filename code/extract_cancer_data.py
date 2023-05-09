import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("projdata.csv", names=["index","inst","time","status","age","sex","ph.ecog","ph.karno","pat.karno","meal.cal","wt.loss"])
print(df)
print(list(df["time"]))
time_col = list(df["time"])
time_col.pop(0)

# only works due to how well-formatted the data is
time_col = [int(time) for time in time_col]
print(time_col)

print(f"max of time is {max(time_col)}")
time_ranges = np.arange(0, 1040, 10)
print(f"{time_ranges=}")

cdf_dict = { time : 0 for time in time_ranges }
print(f"{cdf_dict=}")
# gather a collection of CDF's
for time in time_ranges:
    for col_entry in time_col: 
        if(col_entry <= time):
            cdf_dict[time]+=1
        else:
            continue
    cdf_dict[time] = cdf_dict[time] / len(time_col)

# fig, axs = plt.subplots( figsize=(9, 3), sharey=True)

fig = plt.figure()
ax = plt.subplot(111)

y = [cdf_dict[value] for value in cdf_dict]
print(f"{y=}")
ax.plot(time_ranges, y)
plt.title(f" Plot of F(t) obtained from data versus t")

# print(f"max of tspace is {max(tspace)}")
plt.xlim([min(time_ranges), max(time_ranges)]) # evaluates to ax.xlim([0, 15]) with example data above
plt.xlabel("t (survival time in days)")
plt.ylabel("F(t)")
# plt.legend(legend_entries,  fancybox=True, shadow=True, ncol=5)
plt.show()
# plt.savefig(f"Fwithλ={l}-ν={v}-ω={w}.png")

fig = plt.figure()
ax = plt.subplot(111)
print(cdf_dict)

# if the CDF has been configures properly, the maximum value after a certain
# range should be 1.0, and the minimum value shoule be 0.0

