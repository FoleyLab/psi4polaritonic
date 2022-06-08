"""
Helper function for CQED-CIS in the coherent state basis

References:
    Equations and algorithms from 
    [Haugland:2020:041043], [DePrince:2021:094112], and [McTague:2021:ChemRxiv] 

"""

__authors__ = ["Jon McTague", "Jonathan Foley"]
__credits__ = ["Jon McTague", "Jonathan Foley"]

__copyright_amp__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__ = "BSD-3-Clause"
__date__ = "2021-01-15"

# ==> Import Psi4, NumPy, & SciPy <==
import psi4
import numpy as np
import scipy.linalg as la
import time
from helper_cqed_rhf import cqed_rhf


def cs_cqed_cis(lambda_vector, omega_val, molecule_string, psi4_options_dict):
    """Computes the QED-RHF energy and density

    Arguments
    ---------
    lambda_vector : 1 x 3 array of floats
        the electric field vector, see e.g. Eq. (1) in [DePrince:2021:094112]
        and (15) in [Haugland:2020:041043]

    omega_val : complex float
        the complex energy associated with the photon, see Eq. (3) in [McTague:2021:ChemRxiv]

    molecule_string : string
        specifies the molecular geometry

    psi4_options_dict : dictionary
        specifies the psi4 options to be used in running requisite psi4 calculations

    Returns
    -------
    cqed_cis_dictionary : dictionary
        Contains important quantities from the cqed_rhf calculation, with keys including:
            'RHF ENERGY' -> result of canonical RHF calculation using psi4 defined by molecule_string and psi4_options_dict
            'CQED-RHF ENERGY' -> result of CQED-RHF calculation, see Eq. (13) of [McTague:2021:ChemRxiv]
            'CQED-CIS ENERGY' -> numpy array of complex floats comprising energy eigenvalues of CQED-CIS Hamiltonian
            'CQED-CIS L VECTORS' -> numpy array of complex floats comprising the left eigenvectors of CQED-CIS Hamiltonian

    Example
    -------
    >>> cqed_cis_dictionary = cs_cqed_cis([0., 0., 1e-2], 0.2-0.001j, '''\nMg\nH 1 1.7\nsymmetry c1\n1 1\n''', psi4_options_dictionary)

    """

    # define geometry using the molecule_string
    mol = psi4.geometry(molecule_string)
    # define options for the calculation
    psi4.set_options(psi4_options_dict)
    # run psi4 to get ordinary scf energy and wavefunction object
    # scf_e, wfn = psi4.energy('scf', return_wfn=True)

    # run cqed_rhf method
    cqed_rhf_dict = cqed_rhf(lambda_vector, molecule_string, psi4_options_dict)

    # grab necessary quantities from cqed_rhf_dict
    scf_e = cqed_rhf_dict["RHF ENERGY"]
    cqed_scf_e = cqed_rhf_dict["CQED-RHF ENERGY"]
    wfn = cqed_rhf_dict["PSI4 WFN"]
    C = cqed_rhf_dict["CQED-RHF C"]
    eps = cqed_rhf_dict["CQED-RHF EPS"]
    cqed_rhf_dipole_moment = cqed_rhf_dict["CQED-RHF DIPOLE MOMENT"]

    # Create instance of MintsHelper class
    mints = psi4.core.MintsHelper(wfn.basisset())

    # Grab data from wavfunction

    # number of doubly occupied orbitals
    ndocc = wfn.nalpha()

    # total number of orbitals
    nmo = wfn.nmo()

    # number of virtual orbitals
    nvirt = nmo - ndocc

    # need to update the Co and Cv core matrix objects so we can
    # utlize psi4s fast integral transformation!

    # collect rhf wfn object as dictionary
    wfn_dict = psi4.core.Wavefunction.to_file(wfn)

    # update wfn_dict with orbitals from CQED-RHF
    wfn_dict["matrix"]["Ca"] = C
    wfn_dict["matrix"]["Cb"] = C
    # update wfn object
    wfn = psi4.core.Wavefunction.from_file(wfn_dict)

    # occupied orbitals as psi4 objects but they correspond to CQED-RHF orbitals
    Co = wfn.Ca_subset("AO", "OCC")

    # virtual orbitals same way
    Cv = wfn.Ca_subset("AO", "VIR")

    # 2 electron integrals in CQED-RHF basis
    ovov = np.asarray(mints.mo_eri(Co, Cv, Co, Cv))

    # build the (oo|vv) integrals:
    oovv = np.asarray(mints.mo_eri(Co, Co, Cv, Cv))

    # strip out occupied orbital energies, eps_o spans 0..ndocc-1
    eps_o = eps[:ndocc]

    # strip out virtual orbital energies, eps_v spans 0..nvirt-1
    eps_v = eps[ndocc:]

    # Extra terms for Pauli-Fierz Hamiltonian
    # nuclear dipole
    mu_nuc_x = mol.nuclear_dipole()[0]
    mu_nuc_y = mol.nuclear_dipole()[1]
    mu_nuc_z = mol.nuclear_dipole()[2]

    # l \cdot \mu_nuc for d_c
    l_dot_mu_nuc = lambda_vector[0] * mu_nuc_x
    l_dot_mu_nuc += lambda_vector[1] * mu_nuc_y
    l_dot_mu_nuc += lambda_vector[2] * mu_nuc_z

    # dipole arrays in AO basis
    mu_ao_x = np.asarray(mints.ao_dipole()[0])
    mu_ao_y = np.asarray(mints.ao_dipole()[1])
    mu_ao_z = np.asarray(mints.ao_dipole()[2])

    # transform dipole array to CQED-RHF basis
    mu_cmo_x = np.dot(C.T, mu_ao_x).dot(C)
    mu_cmo_y = np.dot(C.T, mu_ao_y).dot(C)
    mu_cmo_z = np.dot(C.T, mu_ao_z).dot(C)

    # \lambda \cdot < \mu >
    # e.g. line 6 of Eq. (18) in [McTague:2021:ChemRxiv]
    l_dot_mu_exp = 0.0
    for i in range(0, 3):
        l_dot_mu_exp += lambda_vector[i] * cqed_rhf_dipole_moment[i]

    # \lambda \cdot \mu_{el}
    # e.g. line 4 Eq. (18) in [McTague:2021:ChemRxiv]
    l_dot_mu_el = lambda_vector[0] * mu_cmo_x
    l_dot_mu_el += lambda_vector[1] * mu_cmo_y
    l_dot_mu_el += lambda_vector[2] * mu_cmo_z

    # dipole constants to add to E_CQED_CIS,
    #  0.5 * (\lambda \cdot \mu_{nuc})** 2
    #      - (\lambda \cdot <\mu> ) ( \lambda \cdot \mu_{nuc})
    # +0.5 * (\lambda \cdot <\mu>) ** 2
    # Eq. (14) of [McTague:2021:ChemRxiv]
    d_c = (
        0.5 * l_dot_mu_nuc ** 2 - l_dot_mu_nuc * l_dot_mu_exp + 0.5 * l_dot_mu_exp ** 2
    )

    # check to see if d_c what we have from CQED-RHF calculation
    assert np.isclose(d_c, cqed_rhf_dict["DIPOLE ENERGY"])

    # build g matrix and its adjoint
    g = np.zeros((1,ndocc * nvirt), dtype=complex)
    g_dag = np.zeros((ndocc * nvirt, 1), dtype=complex)
    for i in range(0, ndocc):
        for a in range(0, nvirt):
            A = a + ndocc
            ia = i * nvirt + a 
            g[0,ia] = (
                -np.sqrt(omega_val) * l_dot_mu_el[i, A]
            )

    # Now compute the adjoint of g
    g_dag = np.conj(g).T

    # get the A + \Delta matrix and the G matrix, since both 
    # involive <S | H | S> terms where |S> represents singly-excited determinants
    A_p_D = np.zeros((ndocc * nvirt, ndocc * nvirt), dtype=complex)
    G = np.zeros((ndocc * nvirt, ndocc * nvirt), dtype=complex)
    Omega = np.zeros((ndocc * nvirt, ndocc * nvirt), dtype=complex)
    for i in range(0, ndocc):
        for a in range(0, nvirt):
            A = a + ndocc
            ia = i * nvirt + a
            
            for j in range(0, ndocc):
                for b in range(0, nvirt):
                    B = b + ndocc
                    jb = j * nvirt + b
                    
                    # ERI contribution to A + \Delta
                    A_p_D[ia, jb] = (2.0 * ovov[i, a, j, b] - oovv[i, j, a, b])
                    
                    # 2-electron dipole contribution to A + \Delta
                    A_p_D[ia, jb] += 2.0 * l_dot_mu_el[i, A] * l_dot_mu_el[j, B]
                    A_p_D[ia, jb] -= l_dot_mu_el[i, j] * l_dot_mu_el[A, B]
                    
                    # bilinear coupling contributions to G
                    # off-diagonal terms (plus occasional diagonal terms)
                    G[ia, jb] += np.sqrt(omega_val / 2) * l_dot_mu_el[i, j] * (a == b)
                    G[ia, jb] -= np.sqrt(omega_val / 2) * l_dot_mu_el[A, B] * (i == j)
                    
                    # diagonal contributions to A_p_D, G, and \Omega matrix
                    if i == j and a == b:
                        # orbital energy contribution to A + \Delta ... this also includes 
                        # the DSE terms that contributed to the CQED-RHF energy 
                        A_p_D[ia, jb] += eps_v[a] 
                        A_p_D[ia, jb] -= eps_o[i] 
                        
                        # dipole constant energy contribution to A + \Delta
                        A_p_D[ia, jb] += d_c 
                        
                        # diagonal \omega term
                        Omega[ia, jb] = omega_val
                        
                        # diagonal terms (WHICH NEEDS MODIFICATION!  THINK THIS WHOLE BLOCK CAN GO SINCE THE ONLY 
                        # SURVIVNIG DIAGONAL TERMS ARE CAPTURED IMMEDIATEDLY ABOVE)   
                        # expectation value term 
                        G[ia, jb] += np.sqrt(omega_val / 2) * l_dot_mu_exp 
                        
                        # electronic term
                        for k in range(0, ndocc):
                            G[ia, jb] -= np.sqrt(omega_val / 2) * l_dot_mu_el[k, k]
    # define the offsets
    R0_offset = 0
    S0_offset = 1
    R1_offset = ndocc * nvirt + 1
    S1_offset = ndocc * nvirt + 2

    # create Hamiltonian for elements H[ias, jbt]
    Htot = np.zeros((ndocc * nvirt * 2 + 2, ndocc * nvirt * 2 + 2), dtype=complex)

    # build the supermatrix
    # g coupling
    Htot[R0_offset:S0_offset, S1_offset:] = g
    Htot[S0_offset:R1_offset, R1_offset:S1_offset] = g_dag
    Htot[R1_offset:S1_offset, S0_offset:R1_offset] = g
    Htot[S1_offset:,          R0_offset:S0_offset] = g_dag

    # A + \Delta 
    Htot[S0_offset:R1_offset, S0_offset:R1_offset] = A_p_D

    # omega
    Htot[R1_offset, R1_offset] = omega_val

    # A + \Delta + \Omega
    Htot[S1_offset:, S1_offset:] = A_p_D + Omega

    # G coupling
    Htot[S1_offset:,S0_offset:R1_offset] = G
    Htot[S0_offset:R1_offset, S1_offset:] = G

    # create Hamiltonian for elements H[ias, jbt]
    H_TDA_RWA = np.zeros((ndocc * nvirt + 1, ndocc * nvirt + 1), dtype=complex)

    # define the TDA offsets
    TDA_S0_offset = 0
    TDA_R1_offset = ndocc * nvirt


    # build the supermatrix
    # g coupling
    H_TDA_RWA[TDA_R1_offset:, TDA_S0_offset:TDA_R1_offset] = g
    H_TDA_RWA[TDA_S0_offset:TDA_R1_offset, TDA_R1_offset:] = g_dag

    # A + \Delta
    H_TDA_RWA[TDA_S0_offset:TDA_R1_offset, TDA_S0_offset:TDA_R1_offset] = A_p_D

    # omega
    H_TDA_RWA[TDA_R1_offset, TDA_R1_offset] = omega_val

    # diagonalize the total QED-CIS matrix and the 
    ECIS, CCIS = np.linalg.eigh(Htot)

    ETDA, CTDA = np.linalg.eigh(H_TDA_RWA)

    cqed_cis_dict = {
        "RHF ENERGY": scf_e,
        "CQED-RHF ENERGY": cqed_scf_e,
        "CQED-CIS ENERGY": ECIS,
        "CQED-TDA-RWA ENERGY": ETDA,
        "CQED-CIS L VECTORS": CCIS,
    }

    return cqed_cis_dict
