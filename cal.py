from collections import OrderedDict
from sympy import *
from sympy.abc import a,b,c,d,i,j,k,p,q,r,x,y,z,A,D,V,I,W,R,L,C

init_printing()

THDN = symbols('THDN')

equations = OrderedDict([
    ('D_cycle_sq', 2*sqrt(a*b/pi)),
    ('D_hydralic', 2*a*b/(a+b)),     
    ('D_ashrae', 1.3*(a*b)**0.625/(a+b)**0.25),
    ('D_oblong', 1.55*(pi*b*b/4+a*b+b*b)**0.625/(pi*b+2*a-2*b)**0.25),
    ('SINAD', 20*log(1/THDN, 10)),
])


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
    

def adjust_spl_port_area(port_area, speaker_area):
    return 10 * log(port_area/speaker_area, 10).evalf()


def adjust_spl_port_dia(port_dia, speaker_dia):
    return 20 * log(port_dia/speaker_dia, 10).evalf()
    

def cal_sinad(thdn_perc):
    """Calculate SINAD from THD+N in percent"""
    return 20 * log(100/thdn_perc, 10).evalf()

def cal_thdn(sinad):
    """Calculate THD+N in percent from SINAD"""
    return 100/(10**(sinad/20))

