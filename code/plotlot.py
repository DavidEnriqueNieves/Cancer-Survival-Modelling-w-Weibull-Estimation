import numpy as np
import matplotlib.pyplot as plt


import math 
# p is phi
# w is omega
# x is t
# l is lambda
F = lambda x, l, v, p, w : 1 - math.exp(-v * math.pow( (x**l)/(p**l - x** l) , w)  )

h = lambda x, l, v, p, w : ( v * w * l * math.pow(p, l)  * (math.pow(x, (w * l) - 1))/(math.pow( p**l - x**l , w + 1) + 0.001) 
*  math.pow(x**l / (p**l - x **l  + 000.1), w)
)
                            # )


p = 4
l = 2
w = 2
v = 2
tspace = np.linspace(0, p, 50)
# plt.plot()

print(F(2,3, 1, 5, 4))
fig, axs = plt.subplots( figsize=(9, 3), sharey=True)


fig = plt.figure()
ax = plt.subplot(111)
phispace = np.linspace(0, 4, 10)

def plotPhiSpace(phispace):
    for phi in phispace:
        # print(Fspace)
        Fspace = [F(t, l,v,phi,w) for t in tspace]
        hspace = [h(t, l,v,phi,w) if t < phi else 0 for t in tspace]
        # plt.plot(tspace, Fspace)
        ax.plot(tspace, hspace)

    legend_entries = [{"φ={:.2f}".format(phi)} for phi in phispace]
    plt.title(f" Plot of F(t) versus t with λ : {l}, ν : {v}, ω : {w}")
    print(f"max of tspace is {max(tspace)}")
    plt.xlim([min(tspace), max(phispace)]) # evaluates to ax.xlim([0, 15]) with example data above
    plt.xlabel("t")
    plt.ylabel("F")
    plt.legend(legend_entries,   shadow=True, nrows=5)
    # plt.show()
    plt.savefig(f"Fwithλ={l}-ν={v}-ω={w}.png")




plotPhiSpace(phispace)