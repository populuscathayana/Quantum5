"""Local harmonic free-energy differences for the active adsorbate coordinates.

No gas-phase translational/rotational terms and no fictional Cu phonons.
Frequencies are unscaled wavenumbers; remove the one unstable TS mode explicitly.
"""
import math

KB_EV_K = 8.617333262145e-5
HC_EV_CM = 1.239841984332e-4
H_EV_S = 4.135667696e-15
HARTREE_EV = 27.211386245988

def harmonic(frequencies_cm, temperature=298.15, transition_state=False):
    assert temperature>0
    freq=[float(f) for f in frequencies_cm]
    assert len(freq)==3,'Frozen Cu model: exactly three active O coordinates'
    neg=[f for f in freq if f<0]
    assert len(neg)==int(transition_state),'Minimum: zero; transition state: one negative eigenvalue'
    pos=[f for f in freq if f>0]
    assert len(pos)==3-int(transition_state) and all(f>0 for f in pos)
    kT=KB_EV_K*temperature
    eps=[HC_EV_CM*f for f in pos]
    zpe=sum(eps)/2
    occupation_energy=sum(e/math.expm1(e/kT) for e in eps)
    free=zpe+kT*sum(math.log(-math.expm1(-e/kT)) for e in eps)
    entropy=(zpe+occupation_energy-free)/temperature
    return {'temperature_K':temperature,'active_real_modes_cm-1':pos,'excluded_imaginary_modes_cm-1':neg,
            'ZPE_eV':zpe,'thermal_excitation_eV':occupation_energy,'Fvib_eV':free,'Svib_eV_K':entropy,
            'definition':'Local adsorbate harmonic free energy at fixed Cu coordinates and fixed cell; no gas translation/rotation or Cu phonon contribution.'}

def differences(energies_eh, modes, temperature=298.15):
    assert set(energies_eh)==set(modes)=={'hcp','ts','fcc'}
    h={key:harmonic(modes[key],temperature,key=='ts') for key in energies_eh}
    reference=float(energies_eh['hcp'])
    e={key:(float(val)-reference)*HARTREE_EV for key,val in energies_eh.items()}
    f={key:e[key]+h[key]['Fvib_eV']-h['hcp']['Fvib_eV'] for key in e}
    return {'reference':'hcp','electronic_relative_eV':e,'local_harmonic_relative_eV':f,'vibrations':h,
            'electronic_forward_barrier_eV':e['ts'],'electronic_reverse_barrier_eV':e['ts']-e['fcc'],
            'local_harmonic_forward_barrier_eV':f['ts'],'local_harmonic_reverse_barrier_eV':f['ts']-f['fcc'],
            'single_hop_quantum_hTST_rate_s-1':(KB_EV_K*temperature/H_EV_S)*math.exp(-f['ts']/(KB_EV_K*temperature))}
