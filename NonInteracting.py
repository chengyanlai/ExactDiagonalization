#!/usr/bin/env python3
# coding=utf-8
import numpy as np
import scipy.linalg as sla
import os
import copy
import h5py

def Hamiltonian(L):
    elem1 = np.ones(L - 1) * -1.0e0
    return np.diag(elem1, -1) + np.diag(elem1, 1)

def Equation_of_Motion(A, H):
    Commutator = -1.0j * (A.dot(H) - H.dot(A))
    return Commutator

def RK4(Ct, t, dt, Ht):
    """
    My f() is Equation of Motion.

        d A = -i[A, H(t)]

    # k1 = dt * f( x[i], t )
    # k2 = dt * f( x[i] + 0.5 * k1, t + 0.5 * dt )
    # k3 = dt * f( x[i] + 0.5 * k2, t + 0.5 * dt )
    # k4 = dt * f( x[i] + k3, t + dt )
    # x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
    """
    """ Construct H(t) = H1 """
#     H1 = Hf(time=t)

    """ k1 = dt * f( x[t], t ) """
    C1 = dt * ( Equation_of_Motion(Ct, Ht) )

    """ Construct H(t+0.5*dt) = H2 """
#     H2 = Hf(time=t + 0.5 * dt, **kwargs)
    """ Construct x[t] + 0.5 * k1 """
    Cw = Ct + 0.50 * C1
    """ k2 = dt * f( x[i] + 0.5 * k1, t[i] + 0.5 * dt ) """
    C2 = dt * ( Equation_of_Motion(Cw, Ht) )

    """
    H3 = H2
    Construct x[t] + 0.5 * k2
    """
    Cw = Ct + 0.50 * C2
    """ k3 = dt * f( x[i] + 0.5 * k2, t[i] + 0.5 * dt ) """
    C3 = dt * ( Equation_of_Motion(Cw, Ht) )

    """ Construct H(t+dt) = H4 """
#     H4 = Hf(time=t + dt, **kwargs)
    """ Construct x[i] + k3 """
    Cw = Ct + C3
    """ k4 = dt * f( x[i] + k3, t + dt ) """
    C4 = dt * ( Equation_of_Motion(Cw, Ht) )

    """ x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0 """
    return Ct + (C1 + 2.0 * (C2 + C3) + C4) / 6.0

##### Part - 1 #####
# Ls = [121, ]
# Ls = [101, ]
# Ls = [81, 21, ]
# Ls = [61, 41, ]
# for L in Ls:
#     N = L - 1
#
#     H = Hamiltonian(L)
#     Eval, Evec = sla.eigh(H)
#     rho = np.float64(N) * np.outer(Evec[:,0], Evec[:,0])
#
#     T0 = 0.0
#     Tf = 200.0
#     dt = 0.01
#     for LOC in [L - 1,]:#[np.int((L-1)/2), L - 1]:
#         fname = "".join(["L", str(L), "-N", str(N), "-LOC", str(LOC), ".h5"])
#         file = os.path.join("/condensate1/GitRepo/ExactDiagonalization/data/Noninteracting", fname)
#         # file = os.path.join("/Volumes/Files/GitRepo/ExactDiagonalization/data/Noninteracting", fname)
#         if os.path.isfile(file):
#             continue
#         f = h5py.File(file, "w")
#
#         Trs = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
#         dset = f.create_dataset("Trs", data=Trs)
#         Vfinals = [-3.0, -5.0, -7.0, -9.0]
#         dset = f.create_dataset("Vfinals", data=Vfinals)
#
#         for Vfs in Vfinals:
#             Vname = "".join(["V", str(Vfs)])
#             Vg = f.create_group(Vname)
#             for Tr in Trs:
#                 Trname = "".join(["Tr", str(Tr)])
#                 Trg = Vg.create_group(Trname)
#                 CtSoft = copy.copy(rho)
#                 Density = [np.diag(CtSoft).real,]
#                 NS = [CtSoft[LOC,LOC].real,]
#                 time = T0
#                 Ht = copy.copy(H)
#                 Ht[LOC,LOC] = 0.0
#                 # Vdyn = np.linspace(0, Vfs, Tr+1)
#                 while time <= Tf:
#                     if Tr - time > dt:
#                         Ht[LOC,LOC] = Vfs * time / Tr
#                     else:
#                         Ht[LOC,LOC] = Vfs
#                     CtSoft = RK4(CtSoft, time, dt, Ht)
#                     Density = np.vstack([Density, np.diag(CtSoft).real])
#                     NS.append(CtSoft[LOC,LOC].real)
#                     time += dt
#                 Trg.create_dataset("Density", data=Density)
#                 Trg.create_dataset("NS", data=NS)
#             if Vfs == Vfinals[0]:
#                 Tlist = np.linspace(T0, Tf+dt, len(NS))
#         dset = f.create_dataset("Tlist", data=Tlist)
#         f.close()

##### Part - 2 #####
L = 101
N = L - 1

# V1 = -1.0
# Tr1 = 10.0# When t=10t_0, V reach V1. dV/dt = 0.1
# V2 = -4.0# This is the amount increase, not final value. Final value is V1+V2.
# Trps = np.array([13.330, 20.0, 40.0, 80.0, 120.0])

# V1 = -1.0
# Tr1 = 20.0# When t=10t_0, V reach V1. dV/dt = 0.1
# V2 = -4.0# This is the amount increase, not final value. Final value is V1+V2.
# Trps = np.array([13.330, 20.0, 40.0, 80.0, 120.0])

# V1 = -2.0
# Tr1 = 20.0# When t=20t_0, V reach V1. dV/dt = 0.1
# V2 = -3.0# This is the amount increase, not final value. Final value is V1+V2.
# Trps = np.array([10.0, 15.0, 30.0, 60.0, 90.0])

# V1 = -2.0
# Tr1 = 40.0
# V2 = -3.0# This is the amount increase, not final value. Final value is V1+V2.
# Trps = np.array([10.0, 15.0, 30.0, 60.0, 90.0])

# V1 = -3.0
# Tr1 = 30.0
# V2 = -2.0# This is the amount increase, not final value. Final value is V1+V2.
# Trps = np.array([6.660, 10.0, 20.0, 40.0, 60.0])

V1 = -3.0
Tr1 = 60.0
V2 = -2.0# This is the amount increase, not final value. Final value is V1+V2.
Trps = np.array([6.660, 10.0, 20.0, 40.0, 60.0])

H = Hamiltonian(L)
Eval, Evec = sla.eigh(H)
rho = np.float64(N) * np.outer(Evec[:,0], Evec[:,0])

T0 = 0.0
Tf = 200.0
dt = 0.01

for LOC in [np.int((L-1)/2), L - 1]:
    fname = "".join(["L", str(L), "-N", str(N), "-LOC", str(LOC), "-V", str(np.int(np.abs(V1))), "-Tr1", str(np.int(Tr1)), ".h5"])
    file = os.path.join("/condensate1/GitRepo/ExactDiagonalization/data/Noninteracting", fname)
    # file = os.path.join("/Volumes/Files/GitRepo/ExactDiagonalization/data/Noninteracting", fname)
    if os.path.isfile(file):
        continue
    f = h5py.File(file, "w")

    dset = f.create_dataset("Tr1", data=Tr1)
    dset = f.create_dataset("V1", data=V1)
    dset = f.create_dataset("V2", data=V2)
    dset = f.create_dataset("Trps", data=Trps)
    for Trp in Trps:
        print(Trp)
        Trname = "".join(["Trp", str(np.int(Trp))])
        Trg = f.create_group(Trname)

        CtSoft = copy.copy(rho)
        Density = [np.diag(CtSoft).real,]
        NS = [CtSoft[LOC,LOC].real,]

        time = T0
        Ht = copy.copy(H)
        Ht[LOC,LOC] = 0.0
        Vdyn = [0.0,]
        while time <= Tf:
            if time < Tr1:
                Vval = V1 * time / Tr1
            elif np.abs(time-Tr1) < 1.0e-5:
                Vval = V1
            elif time < (Trp + Tr1):
                Vval = V1 + V2 * (time - Tr1) / Trp
            else:
                Vval = V1 + V2
            Vdyn.append(Vval)
            Ht[LOC,LOC] = Vval

            CtSoft = RK4(CtSoft, time, dt, Ht)
            Density = np.vstack([Density, np.diag(CtSoft).real])
            NS.append(CtSoft[LOC,LOC].real)
            time += dt
        Trg.create_dataset("Density", data=Density)
        Trg.create_dataset("NS", data=NS)
        Trg.create_dataset("Vdyn", data=Vdyn)
        if Trp == Trps[0]:
            Tlist = np.linspace(T0, Tf+dt, len(NS))
    dset = f.create_dataset("Tlist", data=Tlist)
    f.close()
