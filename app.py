import math
from statistics import NormalDist

import streamlit as st

st.set_page_config(page_title="Pavement Design", page_icon="🛣️", layout="centered")

st.title("Pavement Design Calculator")
st.caption("Choose a design method, enter your values, and calculate the pavement design.")

METHODS = ["IRC:37-2018 (Flexible Pavement)", "AASHTO 1993"]

def reliability_constants(reliability):
    """Return approximate normal-distribution constants for the IRC fatigue/rutting equations."""
    # The values used in the supplied program were given for 80% and 90%.
    if reliability >= 90:
        return 2.0, 0.85
    return 2.5, 0.85

def calculate_irc(cvpd, cbr, vdf, growth, years, direction, lane_factor,
                  reliability, va, vbe, bituminous_mr, terrain):
    if cvpd <= 0 or cbr <= 0 or vdf <= 0:
        raise ValueError("CVPD, CBR, and VDF must be greater than zero.")
    if years <= 0 or growth < 0:
        raise ValueError("Design period must be positive and growth rate cannot be negative.")
    if not 0 < direction <= 1 or not 0 < lane_factor <= 1:
        raise ValueError("Direction factor and lane distribution factor must be between 0 and 1.")
    if not 0 < reliability < 100:
        raise ValueError("Reliability must be between 0 and 100.")
    if va <= 0 or vbe <= 0 or bituminous_mr <= 0:
        raise ValueError("Va, Vbe, and bituminous resilient modulus must be positive.")

    r = growth / 100.0
    if r == 0:
        cumulative_growth = years
    else:
        cumulative_growth = ((1 + r) ** years - 1) / r

    msa = cvpd * 365 * cumulative_growth * vdf * direction * lane_factor / 1_000_000
    subgrade_mr = 10.0 * (cbr ** 0.64)  # MPa, empirical approximation used in the supplied code

    # Simplified screening approximation; this is not a layered-elastic analysis.
    subgrade_strain = 0.0004 * (msa ** 0.1) / max(subgrade_mr, 1.0) ** 0.25
    bituminous_strain = 0.0002 * (msa ** 0.08) / max(bituminous_mr, 1.0) ** 0.2

    k1, k2 = reliability_constants(reliability)
    fatigue_life = k1 * 1e-4 * (1 / max(bituminous_strain, 1e-9)) ** 3.89 * (1 / max(va, 0.1)) ** 0.854 * (1 / max(vbe, 0.1)) ** 2.63
    rutting_life = k2 * 1e-8 * (1 / max(subgrade_strain, 1e-9)) ** 4.5337

    # Initial thickness estimates and iterative increases follow the structure of the supplied program.
    bc = 40.0
    dbm = 100.0
    wmm = 250.0
    gsb = 200.0
    iterations = 0
    while iterations < 500:
        total_bituminous = bc + dbm
        # Approximate strain adjustment as the bituminous layer gets thicker.
        adjusted_bituminous_strain = bituminous_strain * (100.0 / max(total_bituminous, 1.0)) ** 0.35
        adjusted_subgrade_strain = subgrade_strain * (300.0 / max(wmm + gsb, 1.0)) ** 0.25
        fatigue_life = k1 * 1e-4 * (1 / max(adjusted_bituminous_strain, 1e-9)) ** 3.89 * (1 / max(va, 0.1)) ** 0.854 * (1 / max(vbe, 0.1)) ** 2.63
        rutting_life = k2 * 1e-8 * (1 / max(adjusted_subgrade_strain, 1e-9)) ** 4.5337

        if fatigue_life >= msa * 1_000_000 and rutting_life >= msa * 1_000_000:
            break
        if fatigue_life < msa * 1_000_000:
            dbm += 10.0
            if dbm >= 250:
                bc += 10.0
                dbm = 100.0
        if rutting_life < msa * 1_000_000:
            wmm += 25.0
            if wmm >= 500:
                gsb += 25.0
                wmm = 250.0
        iterations += 1

    return {
        "msa": msa,
        "subgrade_mr": subgrade_mr,
        "subgrade_strain": subgrade_strain,
        "bituminous_strain": bituminous_strain,
        "fatigue_life": fatigue_life,
        "rutting_life": rutting_life,
        "bc": bc,
        "dbm": dbm,
        "wmm": wmm,
        "gsb": gsb,
        "total_thickness": bc + dbm + wmm + gsb,
        "iterations": iterations,
        "terrain": terrain,
    }

def calculate_aashto(cvpd, vdf, reliability, years, growth, direction, lane_factor,
                     cbr, pi, pt, so, a1, a2, a3, m2, m3):
    if cvpd <= 0 or vdf <= 0 or cbr <= 0:
        raise ValueError("CVPD, VDF, and CBR must be greater than zero.")
    if years <= 0 or growth < 0:
        raise ValueError("Design period must be positive and growth rate cannot be negative.")
    if not 0 < direction <= 1 or not 0 < lane_factor <= 1:
        raise ValueError("Direction factor and lane distribution factor must be between 0 and 1.")
    if not 50 < reliability < 100:
        raise ValueError("Reliability must be greater than 50% and less than 100%.")
    if not (pi > pt > 0):
        raise ValueError("Initial serviceability (Pi) must be greater than terminal serviceability (Pt), and Pt must be positive.")

    r = growth / 100.0
    if r == 0:
        growth_factor = years
    else:
        growth_factor = ((1 + r) ** years - 1) / r
    esal = cvpd * 365 * growth_factor * vdf * direction * lane_factor
    mr_psi = 1500.0 * cbr
    zr = NormalDist().inv_cdf(1 - reliability / 100.0)
    delta_psi = pi - pt

    def log10_w18(sn):
        term1 = zr * so
        term2 = 9.36 * math.log10(sn + 1.0) - 0.20
        denominator = 0.40 + (1094.0 / ((sn + 1.0) ** 5.19))
        term3 = math.log10(delta_psi / denominator)
        term4 = 2.32 * math.log10(mr_psi) - 8.07
        return term1 + term2 + term3 + term4

    required_sn = None
    for i in range(1, 2001):
        sn = i * 0.01
        if log10_w18(sn) >= math.log10(max(esal, 1.0)):
            required_sn = sn
            break
    if required_sn is None:
        raise ValueError("Could not find a structural number. Check the input values.")

    # Minimum total layer thickness search, in 0.5-inch steps, preserving the supplied coefficients.
    best = None
    max_inches = 100
    for d1_steps in range(2, max_inches * 2 + 1):
        d1 = d1_steps * 0.5
        for d2_steps in range(2, max_inches * 2 + 1):
            d2 = d2_steps * 0.5
            for d3_steps in range(2, max_inches * 2 + 1):
                d3 = d3_steps * 0.5
                sn = a1 * d1 + a2 * d2 * m2 + a3 * d3 * m3
                if sn >= required_sn:
                    total = d1 + d2 + d3
                    if best is None or total < best[0]:
                        best = (total, d1, d2, d3, sn)
                    break
            if best is not None and d1 + d2 + 1.0 >= best[0]:
                break
        if best is not None and d1 + 1.0 >= best[0]:
            break

    if best is None:
        # A fast practical fallback in case unusual coefficients make the exhaustive search too large.
        d1 = max(2.0, math.ceil((required_sn / max(a1, 0.01)) * 2) / 2)
        remaining = max(0.0, required_sn - a1 * d1)
        d2 = 2.0
        d3 = max(2.0, math.ceil((remaining / max(a3 * m3, 0.01)) * 2) / 2)
        best = (d1 + d2 + d3, d1, d2, d3, a1*d1 + a2*d2*m2 + a3*d3*m3)

    return {
        "esal": esal,
        "mr_psi": mr_psi,
        "zr": zr,
        "required_sn": required_sn,
        "d1": best[1],
        "d2": best[2],
        "d3": best[3],
        "provided_sn": best[4],
        "total_inches": best[0],
        "total_mm": best[0] * 25.4,
    }

method = st.radio("Design method", METHODS, horizontal=False)

if method == METHODS[0]:
    st.subheader("IRC:37-2018 inputs")
    with st.form("irc_form"):
        col1, col2 = st.columns(2)
        with col1:
            cvpd = st.number_input("Commercial vehicles per day (CVPD)", min_value=1, value=1500, step=100)
            cbr = st.number_input("Subgrade CBR (%)", min_value=0.1, value=8.0, step=0.5)
            vdf = st.number_input("Vehicle Damage Factor (VDF)", min_value=0.1, value=3.5, step=0.1)
            growth = st.number_input("Annual traffic growth (%)", min_value=0.0, value=5.0, step=0.5)
            years = st.number_input("Design period (years)", min_value=1, value=20, step=1)
            direction = st.number_input("Directional factor", min_value=0.01, max_value=1.0, value=0.5, step=0.05)
        with col2:
            lane_factor = st.number_input("Lane distribution factor", min_value=0.01, max_value=1.0, value=0.75, step=0.05)
            reliability = st.selectbox("Reliability (%)", [80, 90], index=1)
            va = st.number_input("Air voids, Va (%)", min_value=0.1, value=3.5, step=0.1)
            vbe = st.number_input("Effective bitumen, Vbe (%)", min_value=0.1, value=11.5, step=0.5)
            bituminous_mr = st.number_input("Bituminous resilient modulus (MPa)", min_value=1.0, value=2000.0, step=100.0)
            terrain = st.selectbox("Terrain", ["Plain", "Rolling", "Hilly"], index=0)
        submitted = st.form_submit_button("Calculate pavement design", type="primary")

    if submitted:
        try:
            result = calculate_irc(cvpd, cbr, vdf, growth, years, direction, lane_factor,
                                   reliability, va, vbe, bituminous_mr, terrain)
            st.subheader("IRC design results")
            st.metric("Design traffic", "{:.2f} MSA".format(result["msa"]))
            st.metric("Estimated subgrade resilient modulus", "{:.2f} MPa".format(result["subgrade_mr"]))
            st.write("**Suggested layer thicknesses**")
            st.table({
                "Layer": ["Bituminous Concrete (BC)", "Dense Bituminous Macadam (DBM)", "Wet Mix Macadam (WMM)", "Granular Sub-Base (GSB)"],
                "Thickness (mm)": [result["bc"], result["dbm"], result["wmm"], result["gsb"]],
            })
            st.metric("Total pavement thickness", "{:.0f} mm".format(result["total_thickness"]))
            with st.expander("Calculation details"):
                st.write("Approximate bituminous strain: {:.6g}".format(result["bituminous_strain"]))
                st.write("Approximate subgrade strain: {:.6g}".format(result["subgrade_strain"]))
                st.write("Estimated fatigue life: {:.3g} standard-axle repetitions".format(result["fatigue_life"]))
                st.write("Estimated rutting life: {:.3g} standard-axle repetitions".format(result["rutting_life"]))
                st.write("Iteration count: {}".format(result["iterations"]))
            st.warning("Important: the supplied IRC approach uses simplified strain approximations. Treat the thicknesses as preliminary estimates, not as a validated IRC-compliant design. Verify traffic calculations, material properties, catalogue requirements, and layered-elastic analysis before practical use.")
        except Exception as error:
            st.error(str(error))

else:
    st.subheader("AASHTO 1993 inputs")
    with st.form("aashto_form"):
        col1, col2 = st.columns(2)
        with col1:
            cvpd = st.number_input("Commercial vehicles per day (CVPD)", min_value=1, value=1500, step=100, key="aa_cvpd")
            vdf = st.number_input("Vehicle Damage Factor (VDF)", min_value=0.1, value=3.5, step=0.1, key="aa_vdf")
            reliability = st.number_input("Reliability (%)", min_value=50.1, max_value=99.9, value=90.0, step=1.0)
            years = st.number_input("Design period (years)", min_value=1, value=20, step=1, key="aa_years")
            growth = st.number_input("Annual traffic growth (%)", min_value=0.0, value=5.0, step=0.5, key="aa_growth")
            direction = st.number_input("Directional factor", min_value=0.01, max_value=1.0, value=0.5, step=0.05, key="aa_direction")
            lane_factor = st.number_input("Lane distribution factor", min_value=0.01, max_value=1.0, value=0.8, step=0.05, key="aa_lane")
            cbr = st.number_input("Subgrade CBR (%)", min_value=0.1, value=8.0, step=0.5, key="aa_cbr")
        with col2:
            pi = st.number_input("Initial serviceability (Pi)", min_value=0.1, value=4.2, step=0.1)
            pt = st.number_input("Terminal serviceability (Pt)", min_value=0.1, value=2.5, step=0.1)
            so = st.number_input("Overall standard deviation (So)", min_value=0.1, value=0.45, step=0.05)
            a1 = st.number_input("Layer coefficient a1", min_value=0.01, value=0.44, step=0.01)
            a2 = st.number_input("Layer coefficient a2", min_value=0.01, value=0.14, step=0.01)
            a3 = st.number_input("Layer coefficient a3", min_value=0.01, value=0.11, step=0.01)
            m2 = st.number_input("Drainage coefficient m2", min_value=0.1, value=1.0, step=0.1)
            m3 = st.number_input("Drainage coefficient m3", min_value=0.1, value=1.0, step=0.1)
        submitted = st.form_submit_button("Calculate pavement design", type="primary")

    if submitted:
        try:
            result = calculate_aashto(cvpd, vdf, reliability, years, growth, direction, lane_factor,
                                      cbr, pi, pt, so, a1, a2, a3, m2, m3)
            st.subheader("AASHTO 1993 design results")
            st.metric("Design traffic", "{:,.0f} ESALs".format(result["esal"]))
            st.metric("Required structural number (SN)", "{:.2f}".format(result["required_sn"]))
            st.write("**Selected layer thicknesses**")
            st.table({
                "Layer": ["Asphalt surface layer (D1)", "Base layer (D2)", "Sub-base layer (D3)"],
                "Thickness (inches)": [result["d1"], result["d2"], result["d3"]],
            })
            st.metric("Total thickness", "{:.1f} inches ({:.0f} mm)".format(result["total_inches"], result["total_mm"]))
            with st.expander("Calculation details"):
                st.write("Subgrade resilient modulus: {:.0f} psi".format(result["mr_psi"]))
                st.write("Reliability normal deviate (ZR): {:.3f}".format(result["zr"]))
                st.write("Structural number provided: {:.2f}".format(result["provided_sn"]))
            st.warning("Important: results depend on the assumptions and coefficients entered. Confirm the traffic definition, AASHTO equation implementation, layer coefficients, drainage coefficients, and local design requirements before practical use.")
        except Exception as error:
            st.error(str(error))

st.divider()
st.caption("This app runs on your computer and opens in your web browser. It does not upload your calculations to a website.")
