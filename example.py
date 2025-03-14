import matplotlib.pyplot as plt
from main import *

### Paramaters ###

nSteps = 3 # number of steps for walk. Sequence for walk is aBBO 45 deg -> 0 deg -> 45 deg ->...                  
alphaSq = 0.08 # intensity / mean photon number of coherent state
r = 0.07 # squeezing parameter
eta = 0.07 # overall efficiency
gamma = 0 # phase shift between H,V due to group delay in aBBO
mm = 0.7 # mode matching (HOM visibility)
n_noise = 5e-6 # dark count prob per pump pulse in each mode
max_photons = 2 # maximum number of photons detected

### Run simulation ###

pn_ideal = computeWalkOutput(nSteps, r, alphaSq, eta, gamma, max_photons, n_noise)
pn_imperfect = computeWalkOutputWithMM(nSteps, r, alphaSq, eta, gamma, mm, max_photons, n_noise)


# Keep only H-photon or V-photon modes (tracing-over non detected modes)

pn_ideal_H, pn_ideal_V = traceOverHV(pn_ideal)
pn_imperfect_H, pn_imperfect_V = traceOverHV(pn_imperfect)


### Plotting specific outcomes ###

# look at 1-photon H subspace
oneFolds_ideal = filterProbDict(pn_ideal_H, num_photons=1) 
oneFolds_imperfect = filterProbDict(pn_imperfect_H, num_photons=1) 

# look at 2-photon H subspace
twoFolds_ideal = filterProbDict(pn_ideal_H, num_photons=2) 
twoFolds_imperfect = filterProbDict(pn_imperfect_H, num_photons=2) 

# plot
fig, ax = plt.subplots(figsize = (12,8))

ax.bar(np.arange(len(twoFolds_ideal))+0.1, twoFolds_ideal.values(), color='tab:blue',width=0.2, label='Perfect mode overlap')

ax.bar(np.arange(len(twoFolds_imperfect))-0.1, twoFolds_imperfect.values(), color='tab:orange',width=0.2, label='Imperfect mode overlap')

ax.set_xticks(range(len(twoFolds_imperfect)))

ax.set_xticklabels(list(twoFolds_imperfect.keys()), rotation=65)

plt.title('H-photons walk output')
plt.ylabel('Probability')
plt.xlabel('Detection outcome (t0,t1,t2,...)')

plt.legend()
plt.tight_layout()
plt.show()