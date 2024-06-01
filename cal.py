from termcolor import colored
from collections import OrderedDict
from unit_converter.converter import convert
from sympy import *
from sympy.abc import a,b,c,d,i,j,k,p,q,r,x,y,z,A,D,V,I,W,R,L,C
from sympy.printing import pprint

init_printing()

f3, fb, fc, fs, Vas, Vb, Qtc, Qts, Qes, Qms = symbols('f3 fb fc fs Vas Vb Qtc Qts Qes Qms')

THDN = symbols('THDN')
SOUND_SPEED = 343
print(f'The speed of sound set to {colored(SOUND_SPEED, "green")} m/s\n'
      f'use {colored("set_sound_speed()", attrs=["bold"])} to change the value.')


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
    Eq(fb, fs * ( Vas/(Vb+1) )**0.5),
    Eq(fb, fs * ( sqrt(Vas) / sqrt(Vas/((Qtc/Qts)**2 - 1)+1) )),
    Eq(fb, fs * sqrt( Vb * ((Qtc/Qts)**2 - 1) / (Vb+1))),
    Eq(f3, fb / 2**0.5 * sqrt((1/Qtc**2 - 2) + sqrt((1/Qtc**2-2)**2 + 4))),
]


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

