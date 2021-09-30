# ==> Import Psi4, NumPy, & SciPy <==
import psi4
import numpy as np

# Set Psi4 & NumPy Memory Options
psi4.set_memory('2 GB')

numpy_memory = 2

# rhf/cc-pVDZ optimized geometry of formaldehyde
molstr = """

0 1
O 0.0000000000 0.0000000000 5.91268220e-01
C 0.0000000000 0.0000000000 -5.90400099e-01
H 0.0000000000 9.32184336e-01 -1.17703144e+00
H 0.0000000000 -9.32184336e-01 -1.17703144e+00
no_reorient
symmetry c1
"""

# options dict
options_dict = {'basis': 'cc-pVDZ',
               'save_jk': True, 
               'scf_type': 'pk',
               'e_convergence' : 1e-10,
               'd_convergence' : 1e-10}

psi4.set_options(options_dict)
mol = psi4.geometry(molstr)


energy, wfn = psi4.energy("scf/cc-pVDZ", molecule=mol, return_wfn=True)
