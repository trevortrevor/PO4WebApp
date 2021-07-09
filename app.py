import streamlit as st
import numpy as np
from pHcalc.pHcalc import Acid, Neutral, System
import altair as alt
import pandas as pd


def Li_pH(M_LiOHH2O, M_Li2CO3):
    Carbonate = Acid(pKa=[6.35, 10.33], charge=0, conc=M_Li2CO3)
    Li_carbonate = Neutral(charge=1, conc=M_Li2CO3*2)
    LiOH = Neutral(charge=1, conc=M_LiOHH2O)
    system = System(Carbonate, Li_carbonate, LiOH)
    system.pHsolve()
    return system

def Tri_Di_pH(M_Tri = 0.01, M_Di = 0.01, M_Sulphite = 0.00):
    Tri = Acid(pKa=[2.15, 6.82, 12.38], charge=0, conc=M_Tri)
    NaOH_Tri = Neutral(charge=1, conc=M_Tri*3)
    Di = Acid(pKa=[2.15, 6.82, 12.38], charge=0, conc=M_Di)
    NaOH_Di = Neutral(charge=1, conc=M_Di*2)
    Na_Sulphite = Acid(pKa=[1.857, 7.172], charge=0, conc=M_Sulphite) # pKa Sulphurous Acid = 1.857 & 7.172
    NaOH_Sulphite = Neutral(charge=1, conc=M_Sulphite*2)
    system = System(Tri, NaOH_Tri, Di, NaOH_Di, Na_Sulphite, NaOH_Sulphite)
    system.pHsolve()
    return system

def get_pH_and_conc(ppm_Tri, ppm_Di):
    RMM_Tri = 163.9407
    RMM_Di = 141.9588

    M_Tri = (ppm_Tri / (1000*RMM_Tri))
    M_Di = (ppm_Di / (1000*RMM_Di))
    system = Tri_Di_pH(M_Tri=M_Tri, M_Di=M_Di)
    #ppm_Phosphate = (M_Tri + M_Di) * RMM_Phosphate * 1000
    return system.pH


st.set_page_config(page_title='Chemistry demos', page_icon = None,
                    layout = 'wide', initial_sidebar_state = 'auto')


st.sidebar.title("Options")

LiOH_or_PO4 = st.sidebar.radio("Select application", ["Lithium hydroxide",
            "Phosphates"])

if LiOH_or_PO4 == "Lithium hydroxide":
    rmm_LiOHH2O = 41.96
    rmm_Li2CO3 = 73.89
    # Top section
    st.title("LiOH with carbonate contamination")
    st.write("")
    st.write("")
    st.sidebar.title("Settings")
    pH_or_mass_volume = st.sidebar.radio("Input pH or volume & mass", ["pH",
                "Volume and Mass"])
    if pH_or_mass_volume == "Volume and Mass":
        mass_units = st.sidebar.radio("Units for mass", ["g", "mg"])
        if mass_units == "g":
            mult=1
        else:
            mult = 0.001

        volume = st.sidebar.number_input("Volume / L", min_value=0.0,
                                        max_value=None, value=1.0, step=1.0)
        init_mass_LiOHH2O = mult*st.sidebar.number_input(
                            f"Initial mass of LiOH monohydrate / {mass_units}",
                            min_value=0.0, max_value=None, value=10.0 ,step=1.0)

    else:
        pH = st.sidebar.number_input("Initial expected pH", min_value=7.01,
                        max_value=None, value=11.0, step=1.0)
        volume = 1.0
        init_mass_LiOHH2O = 10**(pH-14) * rmm_LiOHH2O

    values = st.sidebar.slider("Mass conversion percentage", min_value=0.0,
                                    max_value=100.0, value=(0.0, 100.0))

    datapoints = st.sidebar.number_input("Number of datapoints", min_value=1,
                        value=20, step=10)
    show_Li = st.sidebar.checkbox("Show lithium concs", value=False)
    percentages = np.linspace(values[0]/100., values[1]/100., datapoints)

    pH_values = []
    concentrations = []
    for p in percentages:
        M_LiOHH2O=((1-p)*init_mass_LiOHH2O)/rmm_LiOHH2O/volume # Conc LiOHH2O
        M_Li2CO3 = (p*init_mass_LiOHH2O / rmm_LiOHH2O) / (2 * volume)
        system = Li_pH(M_LiOHH2O, M_Li2CO3)
        pH_values.append(system.pH)

        temp = []
        for s in system.species:
            # Check to see if it's the carbonate ion with multiple valencies
            if len(s.alpha(system.pH)) > 1:
                # Get the individual concentrations by multiplying overall conc
                # by the fractional amount of each
                for x in s.conc * s.alpha(system.pH):
                    temp.append(x)
            else:
                # No equilibrium, so just use the overall conc
                temp.append(s.conc)
        concentrations.append(temp)

    pH_values = np.array(pH_values)

    col1, col2 = st.beta_columns(2)

    with col1:
        df = pd.DataFrame(np.vstack([100*percentages, pH_values]).T,
                        columns=["Percentage Conversion", "pH"])

        c = alt.Chart(df).mark_line().encode(
                x='Percentage Conversion', y=alt.Y('pH',
                scale=alt.Scale(zero=False), ),
                tooltip=['Percentage Conversion','pH']).interactive()

        st.altair_chart(c, use_container_width=True)

    with col2:
        concentrations = np.array(concentrations)
        Li_total = (concentrations[:,-2] + concentrations[:,-1]).reshape(-1,1)
        concentrations = np.hstack([100*percentages.reshape(-1, 1),
                                    concentrations, Li_total])
        if not show_Li:
            concentrations = concentrations[:,:4]
            columns = ["Percentage Conversion", "Carbonic acid",
            "Bicarbonate", "Carbonate"]
        else:
            columns = ["Percentage Conversion", "Carbonic acid",
                "Bicarbonate", "Carbonate", "Li from carbonate",
                "Li from hydroxide", "Li Total"]
        df2 = pd.DataFrame(concentrations, columns=columns)
        df2 = df2.melt("Percentage Conversion", var_name="Ion",
                        value_name="Concentration")

        c2 = alt.Chart(df2).mark_line().encode(
            x='Percentage Conversion', y='Concentration', color='Ion',
            tooltip=['Percentage Conversion', 'Ion']).interactive()

        st.altair_chart(c2, use_container_width=True)


else:
    # Top section
    st.title("Phosphates and pH")

    st.sidebar.title("Settings")

    min_ppm_phosphate = st.sidebar.number_input(
                        "Min ppm phosphate",
                        min_value=0.01, max_value=None, value=1.0 ,step=1.0)
    max_ppm_phosphate = st.sidebar.number_input(
                        "Max ppm phosphate",
                        min_value=min_ppm_phosphate, max_value=None,
                        value=100.0, step=1.0)
    show_mix = st.checkbox("Show mixed phosphates", value=False)
    if show_mix:
        mix_ratio = st.slider("Na:PO4 ratio",min_value=2.0, max_value=3.0,
                                    value=2.5)
        mix_ratio = 3.0 - mix_ratio
    n_samples = st.sidebar.number_input(
                        "datapoints",
                        min_value=10, max_value=None, value=25,
                        step=5)
    logspace = st.sidebar.checkbox("Use logarithmic ppm spacing", value=True)
    tri = []
    di = []
    mix = []

    if logspace:
        ppms = np.logspace(np.log10(min_ppm_phosphate),
                        np.log10(max_ppm_phosphate),
                        n_samples)
    else:
        ppms = np.linspace(min_ppm_phosphate,max_ppm_phosphate,n_samples)
    for ppm in ppms:
        tri.append(get_pH_and_conc(ppm, 0))
        di.append(get_pH_and_conc(0, ppm))
        if show_mix:
            mix.append(get_pH_and_conc((1.0 - mix_ratio)*ppm, mix_ratio*ppm))

    tri = np.array(tri)
    di = np.array(di)
    mix = np.array(mix)

    if show_mix:
        df = pd.DataFrame(np.vstack([ppms, tri, di, mix]).T,
                    columns=["ppm phosphate", "Tri", "Di", "Mix"])
    else:
        df = pd.DataFrame(np.vstack([ppms, tri, di]).T,
                    columns=["ppm phosphate", "Tri", "Di"])

    df = df.melt("ppm phosphate", var_name="Reagent", value_name="pH")

    c = alt.Chart(df).mark_line(point=True).encode(
            x='ppm phosphate',
            y=alt.Y('pH:Q', scale=alt.Scale(zero=False)),
            color='Reagent',
            tooltip=['ppm phosphate', 'Reagent', 'pH']
            ).properties(height=600).interactive()

    st.altair_chart(c, use_container_width=True)
