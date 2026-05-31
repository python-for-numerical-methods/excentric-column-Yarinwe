import numpy as np
from scipy import optimize


def column_stress_error(P, L, E, A, r, c, e, sigma_allow):
    angle = (L / (2 * r)) * np.sqrt(P / (E * A))
    secant = 1 / np.cos(angle)

    sigma_max = (P / A) * (1 + (e * c / r**2) * secant)

    return sigma_max - sigma_allow


def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """
    P_critical = optimize.newton(
        lambda P: column_stress_error(P, L, E, A, r, c, e, sigma_allow),
        500000
    )

    return float(P_critical)
