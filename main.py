from numpy import pi
import numpy as np

import strawberryfields as sf
from strawberryfields.ops import S2gate, LossChannel, Coherent, BSgate, ThermalLossChannel


def generate_photon_outcomes(N, max_photons=2):
    
    '''Helper function to create a label for the possible
    detection probabilities of max_photons across N modes.
    
    For example, if N=3 and max_photons=2, we expect outcomes to be:
    {(0,0,0), (0,0,1), (0,1,0), (1,0,0), (0,1,1), (1,1,0), (1,0,1), (2,0,0),
    (0,2,0), (0,0,2)}
    
    Input
    
    int N:             Number of detection modes
    int max_photons:   Maximum number of photons detectable across the N modes
    
    Output
    
    list outcomes:     List of N-tuples which are the detection labels.
    
    '''
    
    outcomes = []

    def generate_outcomes_helper(current_outcome, remaining_photons, current_mode):
        if current_mode == N:
            if remaining_photons >= 0:
                outcomes.append(tuple(current_outcome))
            return None

        for photons_in_mode in range(min(max_photons + 1, remaining_photons + 1)):
            new_outcome = current_outcome + [photons_in_mode]
            generate_outcomes_helper(new_outcome, remaining_photons - photons_in_mode, current_mode + 1)

    generate_outcomes_helper([], max_photons, 0)
    return outcomes


def convolve_probabilities(prob1, prob2, max_photons=2):
    
    '''Helper function to convolve the prob distributions described by 
    prob1 and prob2.
    
    Dictionaries have the form: prob1 = { (1,1,1): 0.23 } where 0.23 is the 
    probability to measure (1,1,1).
    
    Input
    
    dict prob1:        First prob distb (e.g. completely mode matched case)
    dict prob2:        Second prob distb (e.g. completely non-mode matched case)
    int max_photons:   Maximum number of photons detectable across the N modes
    
    
    Output
    
    dict convolved_prob:   Convolution of prob1 and prob2 containing events up
                           to max_photons

    '''
    
    convolved_prob = {}

    for outcome1, prob1_value in prob1.items():
        for outcome2, prob2_value in prob2.items():
            
            # Detection label is given by sum of two labels
            convolved_outcome = tuple(sum(x) for x in zip(outcome1, outcome2))
            
            # If total number of photons detected is within max_photons
            if sum(convolved_outcome) <= max_photons:
                
                # add current value of convolved prob to product of prob1 and
                # prob 2 
                convolved_prob[convolved_outcome] = convolved_prob.get(convolved_outcome, 0) + prob1_value * prob2_value

    return convolved_prob

def normalizeProbDict(p):
    
    '''Helper function which normalizes the prob dictionary p'''
    
    p = {key: value /  sum(p.values()) for key, value in p.items()}
    
    return p
    

def filterProbDict(pDict, num_photons=2):
    
    '''Helper function which returns a normalized probDict only containing the
    num_photons subspace. Useful mainly for plotting.'''
    
    filtered_probabilities = {key: value for key, value in pDict.items() if sum(key) == num_photons}
    
    return normalizeProbDict(filtered_probabilities)

def traceOverHV(pDict):
    
    '''Helper function which returns two probDicts only containing the
    H and V photons subspace by tracing over non-detected photons.'''
    
    H_probabilities = {}
    V_probabilities = {}
    
    for key, value in pDict.items():
        H_outcome = key[::2] # keep every second element starting from idx = 0
        V_outcome = key[1::2] # keep every second element starting from idx = 1
        
        H_probabilities[H_outcome] = H_probabilities.get(H_outcome,0) + value
        V_probabilities[V_outcome] = V_probabilities.get(V_outcome,0) + value
    
    return normalizeProbDict(H_probabilities) , normalizeProbDict(V_probabilities)
    

def computeWalkOutput(nSteps, r, alphaSq, eta, gamma, max_photons, n_noise, etaFock=1):
    
    '''Main function which computes the walk output photon statistics, 
    including most experimental imperfections. Uses strawberryfields in the 
    Gaussian backend.
    
    The output photon statistics are stored in a dictionary pn. 
    An example entry is pn = {(2,0,0,0) : 0.23} where the key (2,0,0,0) is 
    the detection label and value 0.23 is the corresponding probability.
    The prob for a desired detection label (a,b,c,d,...) can be obtained by 
    calling pn[ (a,b,c,d,...) ].
    
    **Detection/mode labels follow this convention:
    
    (a,b,c,d,...) = (H;t0, V;t0, H;t1, V;t1, ...)
    
    The walk circuit can be modified as needed. Note that mode 0 is  
    reserved for the herald in the strwaberyfields circuit. But the herald is 
    not included in the output detection labels. Currently the walk is assumed
    to be aBBO[45deg] -> aBBO[0deg] -> aBBO[45deg] ... where 45 deg is defined
    wrt {H,V} basis.
    
    Input
    
    int nSteps:      Total number of steps for the walk.
    float r:         Squeezing parameter for TMSV source
    float alphaSq:   Coherent state intensity / mean photon number
    float eta:       Efficiency of setup (applies equal loss before all detectors)
    float gamma:     Phase between H and V pol due to BBO
    int max_photons: Max number of photons detected at output of walk.
    float n_noise:   Dark count prob per pump pulse in each mode.
    float eta_fock:  Transmission of heralded photon. Useful for computing walk
                     with imperfect mode matching.
    
    Output
    
    dict pn:         Normalized prob distb containing output walk statistics.
    
    
    '''        
    
    nModes = 2*nSteps + 2 # Start with {|H;t0>,|V;t0>}. 
                          # Each subsequent step introduces 2 new modes
                          
    alpha = np.sqrt(alphaSq)
    
    sf.hbar = 2
    prog = sf.Program(nModes+1)  # the + 1 is due to herald mode
    eng = sf.Engine("gaussian")
    
    with prog.context as q:
        
        # Initializing states input to walk
        # Let 0 be the herald mode
        
        S2gate(r, 0)            | (q[0], q[1])
        LossChannel(etaFock)    | q[1] 
        # etaFock=0 when simulating non-mode-matched prob distb.
        # Otherwise, etaFock=1.
        Coherent(alpha)         | q[2]
        
        
        # Quantum walk 
        
        for stepNumber in range(nSteps+1): # stepNumber 0 does nothing...
            
            # if odd step number, assume aBBO at 45 deg in {H,V} basis
            if stepNumber % 2 == 1: 
                
                for k in range(stepNumber-1, -1, -1): 
                    # Mix modes {H,V} with same time bin (basis change)
                    BSgate(theta=-pi/4, phi=0)  | (q[2*k+1], q[2*k+2])
                    BSgate(theta=-pi/4, phi=0)  | (q[2*k+3], q[2*k+4])
                    
                    # Apply time shift to {V} modes
                    BSgate(theta=pi/2, phi=gamma) | (q[2*k+2], q[2*k+4])
                    
                    # Undo basis change
                    BSgate(theta=pi/4, phi=0)  | (q[2*k+1], q[2*k+2])
                    BSgate(theta=pi/4, phi=0)  | (q[2*k+3], q[2*k+4])
            
            # if even step number, assume aBBO at 0 deg in {H,V} basis
            if stepNumber % 2 == 0:
                
                for k in range(stepNumber-1, -1, -1):       
                    # Apply time shift to {V} modes
                    BSgate(theta=pi/2, phi=gamma) | (q[2*k+2], q[2*k+4])
        
        # mimic a HWP operation when ending walk on aBBO 45 deg                              
        if nSteps % 2 == 1: 
            
            for k in range(nSteps+1): 
                # Mix modes {H,V} with same time bin (basis change)
                BSgate(theta=-pi/4, phi=0)  | (q[2*k+1], q[2*k+2])            
        
           
        # Apply loss + dark counts to all channels (including herald!)        
        for i in range(nModes+1):

            ThermalLossChannel(eta, n_noise)     | q[i]

    
    # Run SF engine
    results = eng.run(prog)
    state = results.state
    
    #mus = state.means()
    #covs = state.cov()        
    
    # Compute vacuum, 1-folds, 2-folds
    
    pn = {}
    
    allLabels = generate_photon_outcomes(nModes, max_photons)
   
    for k in range(len(allLabels)):
        detLabel = [1] + list(allLabels[k]) # add herald in mode 0 for computing pn
        
        # We are assuming projection onto Fock states rather than
        # click detectors... Should be an okay approx as long as <n> << 1.
        
        pn[allLabels[k]] = state.fock_prob(detLabel)
    
    return pn  




def computeWalkOutputWithMM(nSteps, r, alphaSq, eta, gamma, mm, max_photons, n_noise):
    
    '''Function which computes the walk output photon statistics when 
    there is imperfect mode matching (i.e. interference visibility) between the
    heralded photon and coherent state.
    
    To model mode mismatch, we convolve the photon statistics from two "orthogonal"
    (i.e. non-interfering) walks happening simultaenously.
    
    In the first walk (pn_para), we compute pn for a 
    heralded photon and coherent state of intensity mm*alpha**2. 
    
    In the second walk (pn_perp), we compute pn for a vacuum input (in place of
    the heralded photon) and coherent state of intensity (1-mm)*alpha**2 
    [note that we still assume we are waiting for a herald click to get the 
    correct relative probs when convolving. But the heralded photon is removed 
    from the calculation by introducing 100% loss in that mode].
    
    We then assume that our detectors cannot "tell" from which of the two walks 
    the photons came from, and thus convolve their photon statistics. 
    
    Input
    
    int nSteps:      Total number of steps for the walk.
    float r:         Squeezing parameter for TMSV source
    float alphaSq:   Coherent state intensity / mean photon number
    float eta:       Efficiency of setup (applies equal loss before all detectors)
    float gamma:     Phase between H and V pol due to BBO
    float mm:        Mode overlap, i.e. HOM visibility between alpha and Fock
    int max_photons: Max number of photons detected at output of walk.
    float n_noise:   Dark count prob per pump pulse in each mode.

    Output
    
    dict pn_final:   Normalized prob distb containing output walk statistics.
    
    
    '''    
    
    
    pn_para = computeWalkOutput(nSteps, r, mm*alphaSq, eta, gamma, max_photons, n_noise, etaFock=1)
    pn_perp = computeWalkOutput(nSteps, r, (1-mm)*alphaSq, eta, gamma, max_photons, n_noise, etaFock=0)
        
    pn_final = convolve_probabilities(pn_para, pn_perp)

    return normalizeProbDict(pn_final)