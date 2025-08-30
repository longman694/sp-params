import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

rho = st.number_input('Air density (ρ, kg/m^3)', value=1.18)
c = st.number_input('Speed of sound in air (c, m/s)', value=344.7)

Qms = st.number_input('Qms', min_value=0.01, value=4.84, step=0.1)
Qes = st.number_input('Qes', min_value=0.01, value=0.47, step=0.1)

Qts = Qms*Qes / (Qms + Qes)

st.latex(r'Qts = (\frac{1}{Qms} + \frac{1}{Qes})^{-1} = ' + f'{Qts:.2}')

Fs = st.number_input('Fs (Hz)', min_value=1.0, max_value=10000.0, value=74.6)

Re = st.number_input('Re (Ω)', min_value=0.0, value=3.7)
Le_i = st.number_input('Le (mH)', min_value=0.0, value=0.22)
Sd_ = st.number_input('Sd (cm^2)', min_value=0.0, value=56.0)
Vas_ = st.number_input('Vas (l)', min_value=0.0, value=3.1)

Le = Le_i / 1000  # H
Vas = Vas_ / 1000  # m^3
Sd = Sd_ / 100 / 100  # m^2

Cms = Vas/(rho*c*c*Sd*Sd)

st.latex(r'C_{ms} = \frac{Vas}{\rho \cdot c^2 \cdot Sd^2} = ' + f'{Cms*1000:.2}' + r'\;mm/N')

Mms = 1/ (2*np.pi*Fs)**2 / (Cms)  # kg  

st.latex(r'M_{ms} = \frac{1}{{2\pi f_s}^2 \cdot C_{ms}} =' + f'{Mms*1000:.3}' + r'\;g')

Rms = np.sqrt(Mms/Cms)/Qms

st.latex(r'R_{ms} = \frac{1}{Q_{ms}}\sqrt{\frac{M_{ms}}{C_{ms}}} = ' + f'{Rms:.3}' + r'\;N \cdot s/m')

Bl = np.sqrt(Re/Qes * np.sqrt(Mms/Cms))

st.latex(r'Bl = \sqrt{\frac{Re}{Qes} \cdot \sqrt{\frac{Mms}{Cms}}} = ' + f'{Bl:.3}' + r'\;N/A')

#Vb_ = st.number_input('Vb (l)', value=4, min_value=0.0)
Vb_ = np.inf
Vb = Vb_ / 1000

l_over_a = np.inf
#l_over_a = st.number_input('port length/port area', value=1, min_value=0.0)

Res = Bl**2/Rms
Les = Bl**2*Cms
Ces = Mms/Bl**2

Leb = Bl**2 / Sd**2 * Vb/ (rho*c**2)
Cev = Sd**2/Bl**2 * rho * l_over_a

# qes = omegas*res*ces
feq = np.logspace(1, np.log10(20_000), 1_000)
omega = 2 * np.pi * feq

n=1
Zvc = Le*(omega*1j)**n

Re_ = Re + Zvc.real  # freq. dependent resistance
Le_ = Zvc.imag/omega  # freq. dependent inductance

Yacoustic = -1j / (Leb * omega - 1/(omega*Cev))  # Ya = 1/ Za
Zmech = (1/Res + 1/(omega*Les*1j) + omega*Ces*1j + Yacoustic)**(-1)
Ztotal = Zmech + Re_ + 1j*omega*Le_

transferfunc = 1j*(omega*Zmech/Ztotal)*Re_*Ces

eta = Sd**2 * rho/c/2/np.pi/Re_/Ces**2/Bl**2
efficiency = eta * abs(transferfunc)**2
power_spl = 112.1+10 * np.log10(efficiency)


df = pd.DataFrame({
    'feq': feq,
    'impedance': abs(Ztotal),
    'phase': np.angle(Ztotal)*180/np.pi,
    'spl': power_spl, 
    'gd': -np.gradient(np.unwrap(np.angle(transferfunc)))/np.gradient(omega)*1000,
})


n0 = Sd**2 * rho/c/2/np.pi/Re/Ces**2/Bl**2

st.latex(r'\eta = \frac{Sd^2}{2\pi\rho c \cdot R_e \cdot Ces^2 \cdot (Bl)^2} = ' + f'{n0 * 100: 0.4} \%')


st.write("---")

st.latex(r'R_{es} = \frac{(Bl)^2}{R_{ms}} =' + f'{Res:.3}' + r'\;\Omega')
st.latex(r'L_{es} = \frac{(Bl)^2}{C_{ms}} =' + f'{Les*1000:.3}' + r'\;mH')
st.latex(r'C_{es} = \frac{M_{ms}}{(Bl)^2} =' + f'{Ces*1000000}' + r'\;uF')



imp_chart = make_subplots(specs=[[{"secondary_y": True}]])
imp_chart.add_trace(
    go.Scatter(x=df['feq'], y=df['impedance'], name="Impedance"),
    secondary_y=False,
)

imp_chart.add_trace(
    go.Scatter(x=df['feq'], y=df['phase'], name="Phase"),
    secondary_y=True,
)

imp_chart.update_xaxes(title_text="Frequency", type="log", range=[1, 4])
magnitude_chart = px.line(df, x='feq', y='spl', log_x=True, range_x=[10, 10000])
group_delay_chart = px.line(df, x='feq', y='gd', log_x=True, range_x=[10, 10000])

tab1, tab2, tab3 = st.tabs(["Impedance", "Magnitude", "Group Delay (ms)"])
with tab1:
    st.plotly_chart(imp_chart, use_container_width=True)
with tab2:
    st.plotly_chart(magnitude_chart, key='magnitude', use_container_width=True)
with tab3:
    st.plotly_chart(group_delay_chart, key='group_delay', use_container_width=True)
