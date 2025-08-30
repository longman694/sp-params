from typing import Literal

import numpy as np

from termcolor import colored
from collections import OrderedDict
from unit_converter.converter import convert
from sympy import *
from sympy.abc import a,b,c,d,i,j,k,p,q,r,x,y,z,A,D,V,I,W,R,L,C
from sympy.printing import pprint

init_printing()

A, V, L = symbols('A V L')
f3, fb, fc, fs, Vas, Vb, Qtc, Qts, Qes, Qms = symbols('f3 fb fc fs Vas Vb Qtc Qts Qes Qms')
k, Dv, Np, Lv = symbols('k Dv Np Lv')

THDN = symbols('THDN')
SOUND_SPEED = 343
AIR_DENSITY = 1.18  # 1.225
print(f'The speed of sound set to {colored(SOUND_SPEED, "green")} m/s\n'
      f'The air density set to {colored(AIR_DENSITY, "green")} kg/m^3\n'
      f'use {colored("set_sound_speed()", attrs=["bold"])} and {colored("set_air_density()", attrs=["bold"])} '
      'to change the value.')


port_eqs = OrderedDict([
    ('D_cycle_sq', Eq(D, 2*sqrt(a*b/pi))),
    ('D_hydralic', Eq(D, 2*a*b/(a+b))),     
    ('D_ashrae', Eq(D, 1.3*(a*b)**0.625/(a+b)**0.25)),
    ('D_oblong', Eq(D, 1.55*(pi*b*b/4+a*b+b*b)**0.625/(pi*b+2*a-2*b)**0.25)),
])

filter_eqs = OrderedDict([
    ('RC 1st order high pass filter', Eq(fc, 1/(2*pi*R*C))),
    ('LC 2nd order high pass filter', Eq(fc, 1/(2*pi*sqrt(L*C)))),
])

sealed_box_eqs = [
    Eq(Vb, Vas / ((Qtc/Qts)**2 - 1)),
    Eq(fb, fs * Qtc / Qts),
    Eq(fb, fs * ( sqrt(Vas) / sqrt(Vas/((Qtc/Qts)**2 - 1)+1) )),
    Eq(fb, fs * sqrt( Vb * ((Qtc/Qts)**2 - 1) / (Vb+1))),
    Eq(f3, fb / 2**0.5 * sqrt((1/Qtc**2 - 2) + sqrt((1/Qtc**2-2)**2 + 4))),
]

ported_box_eqs = OrderedDict([
    ('Fb (Dv)', Eq(fb, c/2/pi * sqrt( pi*Dv*Dv/4*Np / Vb / (Lv + k*Dv) )))
    ('Fb (A)', Eq(fb, c/2/pi * sqrt(A*Np / Vb / (Lv + k*(sqrt(4/pi*A))) )))
    ('Lv (Dv)', Eq(10*c*c/16/pi * Dv*Dv*Np/Vb/fb/fb - k*Dv))
])

k_values = OrderedDict([
    ('free-free', 0.614),
    ('free-flanged', 0.732),
    ('flanged-flanged', 0.850),
    ('bottom slot port', 2.227),
])

helmholtz_resonator_eq = Eq(fb, SOUND_SPEED/2/pi * sqrt(A/V/L))

def print_eqs(eqs):
    if isinstance(eqs, OrderedDict):
        for name, eq in eqs.items():
            print(f'----[{name}]-----')
            pprint(eq)
            print('\n\n')
        return
    for i, eq in enumerate(eqs):
        print(f'-----[{i+1}]-----')
        pprint(eq)
        print('\n\n')


def set_sound_speed(speed):
    global SOUND_SPEED
    SOUND_SPEED = speed
    print(f"Set speed of sound to {SOUND_SPEED} m/s")


def set_air_density(d):
    global AIR_DENSITY
    AIR_DENSITY = d
    print(f"Set air density to {AIR_DENSITY} kg/m^3")


def cal_dim(vol):
    ratios = [
        (1.0, 1.17, 1.47, ''),
        (1.0, 1.45, 2.1, 'tallest'),
        (1.0, 1.28, 1.54, ''),
        (1.0, 1.26, 1.59, 'golden ratio'),
        (1.0, 1.14, 1.39, 'most cubic'),
        (1.0, 1.12, 1.41, 'bed fittest'),
        (1.0, 1.6, 2.6, 'LDC'),
        (1.0, 1.144, 2.0, 'LDC2'),
    ]
    product = []
    for r in ratios:
        product = r[0]*r[1]*r[2]
        unit = (vol/product) ** (1/3)
        print('{:13s}: {:0.2f} {:0.2f} {:0.2f}'.format(r[3], r[0]*unit, r[1]*unit, r[2]*unit))
    


def format_length(meter):
    if meter >= 1000:
        return f'{meter/1000:.3f} km'
    elif 1 > meter >= 0.01:
        return f'{meter*100:.3f} cm'
    elif meter < 0.01:
        return f'{meter*1000:.3f} mm'
    return f'{meter:.3f} m'


def cal_wave_length(frequency):
    """return wave length in meter"""
    wave_length = SOUND_SPEED / frequency
    print(f'  λ = {format_length(wave_length)}')
    print(f'λ/2 = {format_length(wave_length/2)}')
    print(f'λ/4 = {format_length(wave_length/4)}')
    return wave_length


def adjust_spl_port_area(port_area, speaker_area):
    """Add this value to the port SPL (distract if this value is negative)"""
    return 10 * log(port_area/speaker_area, 10).evalf()


def adjust_spl_port_dia(port_dia, speaker_dia):
    """Add this value to the port SPL (distract if this value is negative)"""
    return 20 * log(port_dia/speaker_dia, 10).evalf()
    

def cal_freq_from_len(length):
    return SOUND_SPEED/length


def cal_freq_from_half_len(half_len):
    return SOUND_SPEED/half_len/2


def cal_freq_from_quarter_len(quarter_len):
    return SOUND_SPEED/quarter_len/4


def cal_sinad(thdn_perc):
    """Calculate SINAD from THD+N in percent"""
    return 20 * log(100/thdn_perc, 10)


def cal_thdn(sinad):
    """Calculate THD+N in percent from SINAD"""
    return 100/(10**(sinad/20))


# def cal_spl_sd_xmax(xmax, sd, freq):
#     """xmax in mm, Sd in cm^2"""
#     sd = sd * 0.0001
#     return 20 * log(xmax * sd * freq * freq * 2 * SOUND_SPEED * AIR_DENSITY / pi, 10).evalf()


def cal_spl_sd_xmax(xmax, sd, freq):
    """xmax in mm, Sd in cm^2"""
    xmax = xmax * 0.001
    sd = sd * 0.0001
    vd = sd * xmax
    return 112 + 10 * log(4 * pi**3 * AIR_DENSITY / SOUND_SPEED * vd**2 * freq**4, 10).evalf()


class Crossover:
    """uncompleted"""
    f0, C1, C2, L1, L2 = symbols('f0 C1 C2 L1 L2')

    # low pass filter
    L_EQ = [
        Eq(f0, 1 / pi / R / L1),
        Eq(f0, 1 / pi / sqrt(C1 * L1)),
        Eq(f0, 1 / pi / pow(C1 * L1 * L2, 0.333333)),
        Eq(f0, 1 / pi / pow(C1 * C2 * L1 * L2, 0.25)),
    ]
    # high pass filter
    H_EQ = [
        Eq(f0, 1 / pi / R / C1),
        Eq(f0, 1 / pi / sqrt(C1 * L1)),
        Eq(f0, 1 / pi / pow(C1 * C2 * L1, 0.333333)),
        Eq(f0, 1 / pi / pow(C1 * C2 * L1 * L2, 0.25)),
    ]

    def __init__(self, side: Literal['high', 'low'], impedance, order=2, q=0.707):
        self.R = impedance
        self.Q = q
        if side.lower() == 'high':
            self.eq = self.H_EQ[order - 1]
        elif side.lower() == 'low':
            self.eq = self.L_EQ[order - 1]
        else:
            raise ValueError(f'side must be either "high" or "low"')

    def cal_L_C(self, freq):
        return

