import ctypes as ct
import numpy as np
from . import ROOT

MASS = np.array([.0219818, .0380903, .038064, .022289, .036461, .09808, .0620271])
A7 = ct.c_double * 7


class Cell(ct.Structure):
    _fields_ = [('aqueous', A7), ('mineral', A7), ('water', ct.c_double)]


class Result(ct.Structure):
    _fields_ = [(k, ct.c_double) for k in ['rho', 'water', 'mass', 'volume', 'pH', 'aw']] + [
        ('si', A7), ('source', A7), ('change', A7), ('mass_source', ct.c_double),
        ('accepted', ct.c_int), ('rejected', ct.c_int), ('evaluations', ct.c_int)]


class Chemistry:
    def __init__(self, temperature_c=25., pressure_atm=1.):
        self.lib = ct.CDLL(str(ROOT/'bin/libSaltChemistry.so'))
        self.handle = ct.c_void_p()
        self.lib.salt_chemistry_create.argtypes = [ct.c_char_p, ct.c_double, ct.c_double, ct.POINTER(ct.c_void_p)]
        self.lib.salt_chemistry_initialize.argtypes = [ct.c_void_p, A7, ct.c_double, A7, ct.POINTER(Cell), ct.POINTER(Result)]
        self.lib.salt_chemistry_destroy.argtypes = [ct.c_void_p]
        self.lib.salt_chemistry_last_error.restype = ct.c_char_p
        self.check(self.lib.salt_chemistry_create(str(ROOT/'src/chemistry/pitzer.dat').encode(), temperature_c+273.15, pressure_atm, ct.byref(self.handle)))

    def check(self, code):
        if code:
            raise RuntimeError(self.lib.salt_chemistry_last_error().decode())

    def fractions(self, totals):
        cell, result = Cell(), Result()
        self.check(self.lib.salt_chemistry_initialize(self.handle, A7(*totals), 1., A7(*([0.]*7)), ct.byref(cell), ct.byref(result)))
        return MASS*np.array(cell.aqueous)/result.mass, result

    def saturated_halite(self):
        lo, hi = .01, 15.
        if self.fractions([lo,0,0,0,lo,0,0])[1].si[0]>=0 or self.fractions([hi,0,0,0,hi,0,0])[1].si[0]<=0:
            raise ValueError('Halite saturation is outside the supported inlet bracket')
        for _ in range(48):
            mid=(lo+hi)/2
            if self.fractions([mid,0,0,0,mid,0,0])[1].si[0]>0:
                hi=mid
            else:
                lo=mid
        molality=(lo+hi)/2
        return self.fractions([molality,0,0,0,molality,0,0])[0]

    def close(self):
        if self.handle:
            self.lib.salt_chemistry_destroy(self.handle)
            self.handle=ct.c_void_p()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
