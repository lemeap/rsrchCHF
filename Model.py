# Import Modules
import time
from CoolProp.CoolProp import PropsSI
from PhysicalProperty import *
from Numeric import *
import sympy as sp


class Model(PhysicalProperty):
    def __init__(self):
        print("Models of boiling heat transfer is successfully started.")

    # Models or correlations
    def cal_SZ(self, q, rhof, dh, g, cpf, kf, Pe, lam):  # Saha and Zuber (1974)
        """
        :param q: Heat flux [MW/m2]
        :param rhof: Liquid density [kg/m3]
        :param dh: Hydraulic diameter [m]
        :param g: Mass flux [kg/m2-s]
        :param cpf: Liquid specific heat [J/kg-K]
        :param kf: Thermal conductivity [W/m-K]
        :param Pe: Peclet number [-]
        :param lam: Heat of vaporization [J/kg]
        :return: Equilibrium thermal quality [-]
        """
        if Pe <= 70000:
            dtOSV = 0.0022 * (q * (10 ** 6) * dh) / kf
            xOSV = -cpf * dtOSV / (lam * 10 **3) 
            return round(dtOSV,6), round(xOSV,6)
        else:
            dtOSV = 153.8 * (q * (10 ** 6)) / (g * cpf)
            xOSV = -cpf * dtOSV / (lam * 10 **3)
            return round(dtOSV,6), round(xOSV,6)

    def cal_PSZ(self, q, rhof, dh, g, cpf, kf, Pe, lam):  # Park (2004)
        """
        :param q: Heat flux [MW/m2]
        :param rhof: Liquid density [kg/m3]
        :param dh: Hydraulic diameter [m]
        :param g: Mass flux [kg/m2-s]
        :param cpf: Liquid specific heat [J/kg-K]
        :param kf: Thermal conductivity [W/m-K]
        :param Pe: Peclet number [-]
        :param lam: Heat of vaporization [J/kg]
        :return: Equilibrium thermal quality [-]
        """
        if Pe <= 70000:
            dtOSV = 0.0022 * (q * (10 ** 6) * dh) / kf
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)
        elif Pe > 200000:
            dtOSV = (
                0.08923 * np.exp(-(Pe * (kf ** 0.45) / (dh ** 0.53 * g ** 0.37)) / 25313.63287)
                + 0.00659 * np.exp(-(Pe * (kf ** 0.45) / (dh ** 0.53 * g ** 0.37)) / 211422.70151)
                + 0.00146
            )
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)
        else:
            dtOSV = 153.8 * (q * (10 ** 6)) / (g * cpf)
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)

    def cal_SZde(
        self, q, rhof, dh, g, cpf, kf, Pe, lam, hsur, geo, doi, dio, lh
    ):  # Modified Saha and Zuber (2013)
        """
        :param q: Heat flux [MW/m2]
        :param rhof: Liquid density [kg/m3]
        :param dh: Hydraulic diameter [m]
        :param g: Mass flux [kg/m2-s]
        :param cpf: Liquid specific heat [J/kg-K]
        :param kf: Thermal conductivity [W/m-K]
        :param Pe: Peclet number [-]
        :param lam: Heat of vaporization [J/kg]
        :param hsur: Heated surface of cross-sectional geometry of channel [-]
        :param geo: Channel geometry [-]
        :param doi: Outer or longer length of channel geometry [m]
        :param dio: Inner or shorter length of channel geometry [m]
        :param lh: Heated length of channel [m]
        :return: Equilibrium thermal quality [-]
        """
        qq = q * (10 ** 6)
        R_phpw_1h = doi / (2 * (doi + dio))
        R_phpw_2h = doi / (doi + dio)
        A_phpw_1h = dio / (2 * (doi + dio))
        A_phpw_2h = doi / (doi + dio)

        if Pe > 70000:
            if geo == "R":
                if hsur == 1:
                    dtOSV = ((153.8 * (qq)) / (g * cpf)) * R_phpw_1h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
                else:
                    dtOSV = ((153.8 * (qq)) / (g * cpf)) * R_phpw_2h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
            elif geo == "A":
                if hsur == 1:
                    dtOSV = ((153.8 * (qq)) / (g * cpf)) * A_phpw_1h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
                else:
                    dtOSV = ((153.8 * (qq)) / (g * cpf)) * A_phpw_2h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
            else:
                dtOSV = 153.8 * (qq) / (g * cpf)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        else:
            if geo == "R":
                if hsur == 1:
                    dtOSV = (0.0022 * (qq * dh) / kf) * R_phpw_1h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
                else:
                    dtOSV = (0.0022 * (qq * dh) / kf) * R_phpw_2h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
            elif geo == "A":
                if hsur == 1:
                    dtOSV = (0.0022 * (qq ** dh) / kf) * A_phpw_1h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
                else:
                    dtOSV = 0.0022 * (qq * dh) / kf
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV,6), round(xOSV,6)
            else:
                dtOSV = 0.0022 * (qq * dh) / kf
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)

    def cal_Levy(self, sigma, dh, rhof, muf, kf, re, pr, cpf, g, q, lam, v):  # Levy (1967)
        """
        :param sigma: Suface tension [N/m]
        :param dh: Hydraulic diameter [m]
        :param rhof: Liquid density [kg/m3]
        :param muf: Liquid viscosity [kg/m-s]
        :param kf: Thermal conductivity [W/m-K]
        :param re: reynolds number [-]
        :param pr: prandtl number [-]
        :param cpf: Liquid specific heat [J/kg-K]
        :param g: Mass flux [kg/m2-s]
        :param q: Heat flux [MW/m2]
        :param lam: Heat of vaporization [J/kg]
        :param v: flow velocity [m/s]
        :return: Equilibrium thermal quality [-]
        """
        YB = 0.015 * np.sqrt(sigma * dh * rhof) / muf
        h = 0.023 * (kf / dh) * (re ** 0.8) * (pr ** 0.4)
        f = 0.0055 * (1 + (20000 * (10 ** -4) + 10 ** 6 / re) ** (1 / 3))
        tau = (f * g ** 2) / (8 * rhof)
        Q = (q * 10 ** 6) / (rhof * cpf * np.sqrt(tau / rhof))
        if YB > 30:
            dtOSV = (q * 10 ** 6) / h - 5 * Q * (pr + np.log(1 + 5 * pr) + 0.5 * np.log(YB / 30))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)
        elif YB <= 5:
            dtOSV = (q * 10 ** 6) / h - Q * pr * YB
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)
        else:
            dtOSV = (q * 10 ** 6) / h - 5 * Q * (pr + np.log(1 + pr * ((YB / 5) - 1)))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)

    def cal_Bowr(self, p, q, v, lam, cpf):  # Bowring (1960)
        """
        :param p: pressure [bar]
        :param q: Heat flux [MW/m2]
        :param v: flow velocity [m/s]
        :param lam: Heat of vaporization [J/kg]
        :param cpf: Liquid specific heat [J/kg-K]
        :return: Equilibrium thermal quality [-]
        """
        dtOSV = ((14 + p / 10)) * q / v
        xOSV = -cpf * dtOSV / lam
        return round(dtOSV,6), round(xOSV,6)

    def cal_Unal(self, q, pr, dh, v, cpf, kf, re, refri, lam):  # Unal (1975)
        """
        :param q: Heat flux [MW/m2]
        :param pr: prandtl number [-]
        :param dh: Hydraulic diameter [m]
        :param v: flow velocity [m/s]
        :param cpf: Liquid specific heat [J/kg-K]
        :param kf: Thermal conductivity [W/m-K]
        :param re: reynolds number [-]
        :param refri: refrigerants
        :param lam: Heat of vaporization [J/kg]
        :return: Equilibrium thermal quality [-]
        """
        if refri == "Water":
            if v > 0.45:
                dtOSV = (0.24 * (q * 10 ** 6)) / ((kf / dh) * 0.023 * re ** 0.8 * pr ** 0.4)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            elif v <= 0.45:
                dtOSV = (0.11 * (q * 10 ** 6)) / ((kf / dh) * 0.023 * re ** 0.8 * pr ** 0.4)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        else:
            if v > 0.45:
                dtOSV = (0.28 * (q * 10 ** 6)) / ((kf / dh) * 0.023 * re ** 0.8 * pr ** 0.4)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            elif v <= 0.45:
                dtOSV = (0.11 * (q * 10 ** 6)) / ((kf / dh) * 0.023 * re ** 0.8 * pr ** 0.4)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)

    def cal_Jeong(self, geo, q, dh, dl, v, rhof, rhov, muf, muv, cpf, kf, lam, sigma, Pe, Pr):  # Jeong and Shim (2021)
        """
        :param q: Heat flux [MW/m2]
        :param rhof: Liquid density [kg/m3]
        :param dh: Hydraulic diameter [m]
        :param v: flow velocity [m/s]
        :param cpf: Liquid specific heat [J/kg-K]
        :param kf: Thermal conductivity [W/m-K]
        :param Pe: Peclet number [-]
        :param lam: Heat of vaporization [J/kg]
        :return: Equilibrium thermal quality [-]
        """
        lc = np.sqrt(sigma/(9.8*(rhof-rhov)))
        Re = rhof*v*dh/muf
        We = np.sqrt(rhof*v**2*dl/sigma)
        if We <= 10:
            #dtOSV = ((v ** 2)/ (2*rhof*( 0.8877* (Ca/Bo) ** (1.155))))*10**6
            #dtOSV = (q * 10 ** 6) / (cpf * rhof * v * (( 190 * Re ** -0.91 * Pr ** -0.91)))
            dtOSV = (q * 10 ** 6) / (cpf * rhof * v * (( 450 / (Pe))*(dl/dh) ** -0.25)) #2023/02/16
            #dtOSV = (q * 10 ** 6) / (cpf * rhof * v * (( 370 / (Re * Pr) ))) #2023/02/20
            #dtOSV = (q * (10 ** 6) * dh) / (kf*(865 * (Bo) ** 0.075))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,4), round(xOSV,6), round(We, 4)
        else:
            #dtOSV = (q * 10 ** 6 ) / (cpf * rhof * v * (0.122 * Re ** -0.25 * Pr ** -0.25))
            dtOSV = (q * 10 ** 6 ) / (cpf * rhof * v * (0.0705 * Re ** -0.20 * Pr ** -0.20)) #2023/02/16
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,4), round(xOSV,6), round(We, 4)

    def calLee(self, q, doi, lh, cpf, dtin, dio, geo, hsur, dh, gsat, lam):
        """

        :param q: Heat flux [MW/m2]
        :param doi: Outer or longer length of channel geometry [m]
        :param lh: Heated length of channel [m]
        :param cpf: Liquid specific heat [J/kg-K]
        :param dtin: Inlet liquid subcooling [K]
        :param dio: Inner or shorter length of channel geometry [m]
        :param geo: Channel geometry [-]
        :param hsur: Heated surface of cross-sectional geometry of channel [-]
        :param dh: Hydraulic diameter [m]
        :param gsat: gsat: Saturation mass flux [kg/m2-s]
        :param lam: Heat of vaporization [J/kg]
        :return: Equilibrium thermal quality [-]
        """
        qq = q * (10 ** 6)
        R_heated_1h = doi * lh
        R_heated_2h = 2 * doi * lh
        R_flow = doi * dio
        A_heated_1h = (np.pi * dio) * lh
        A_heated_2h = (np.pi * (dio + doi)) * lh
        A_heated_3h = np.pi * doi
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        C_heated = (np.pi) * doi * lh
        C_flow = (np.pi / 4) * doi ** 2
        R_phpw_1h = doi / (2 * (doi + dio))
        R_phpw_2h = doi / (doi + dio)
        A_phpw_1h = dio / (2 * (doi + dio))
        A_phpw_3h = doi / (doi + dio)

        if geo == "R":
            if hsur == 1:
                dtOSV = (qq * R_heated_1h) / (R_flow * cpf * ((gsat + 27) / 0.58))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            else:
                dtOSV = (qq * R_heated_2h) / (R_flow * cpf * ((gsat + 27) / 0.58))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        elif geo == "A":
            if hsur == 1:
                dtOSV = (qq * A_heated_1h) / (A_flow * cpf * ((gsat + 27) / 0.58))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            else:
                dtOSV = (qq * A_heated_2h) / (A_flow * cpf * ((gsat + 27) / 0.58))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        else:
            dtOSV = (qq * C_heated) / (C_flow * cpf * ((gsat + 27) / 0.58))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)

    def calKennedy(self, q, doi, lh, cpf, dtin, dio, geo, hsur, dh, gsat, lam):
        """
        :param q: Heat flux [MW/m2]
        :param doi: Outer or longer length of channel geometry [m]
        :param lh: Heated length of channel [m]
        :param cpf: Liquid specific heat [J/kg-K]
        :param dtin: Inlet liquid subcooling [K]
        :param dio: Inner or shorter length of channel geometry [m]
        :param geo: Channel geometry [-]
        :param hsur: Heated surface of cross-sectional geometry of channel [-]
        :param dh: Hydraulic diameter [m]
        :param gsat: gsat: Saturation mass flux [kg/m2-s]
        :param lam: Heat of vaporization [J/kg]
        :return: Equilibrium thermal quality [-]
        """
        qq = q * (10 ** 6)
        R_heated_1h = doi * lh
        R_heated_2h = 2 * doi * lh
        R_flow = doi * dio
        A_heated_1h = (np.pi * dio) * lh
        A_heated_2h = (np.pi * (dio + doi)) * lh
        A_heated_3h = np.pi * doi
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        C_heated = (np.pi) * doi * lh
        C_flow = (np.pi / 4) * doi ** 2
        R_phpw_1h = doi / (2 * (doi + dio))
        R_phpw_2h = doi / (doi + dio)
        A_phpw_1h = dio / (2 * (doi + dio))
        A_phpw_3h = doi / (doi + dio)

        if geo == "R":
            if hsur == 1:
                dtOSV = (qq * R_heated_1h) / (R_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            else:
                dtOSV = (qq * R_heated_2h) / (R_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        elif geo == "A":
            if hsur == 1:
                dtOSV = (qq * A_heated_1h) / (A_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            elif hsur == 2:
                dtOSV = (qq * A_heated_2h) / (A_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            elif hsur == 3:
                dtOSV = (qq * A_heated_2h) / (A_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        else:
            dtOSV = (qq * C_heated) / (C_flow * cpf * (gsat * 1.11))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)

    def calAl_Yahia(self, q, doi, lh, cpf, dtin, dio, geo, hsur, dh, gsat, lam, p):
        """
        :param q: Heat flux [MW/m2]
        :param doi: Outer or longer length of channel geometry [m]
        :param lh: Heated length of channel [m]
        :param cpf: Liquid specific heat [J/kg-K]
        :param dtin: Inlet liquid subcooling [K]
        :param dio: Inner or shorter length of channel geometry [m]
        :param geo: Channel geometry [-]
        :param hsur: Heated surface of cross-sectional geometry of channel [-]
        :param dh: Hydraulic diameter [m]
        :param gsat: Saturation mass flux [kg/m2-s]
        :param lam: Heat of vaporization [J/kg]
        :param p: pressure [bar]
        :return: Equilibrium thermal quality [-]
        """
        qq = q * (10 ** 6)
        R_heated_1h = doi * lh
        R_heated_2h = 2 * doi * lh
        R_flow = doi * dio
        A_heated_1h = (np.pi * dio) * lh
        A_heated_2h = (np.pi * (dio + doi)) * lh
        A_heated_3h = np.pi * doi
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        C_heated = (np.pi) * doi * lh
        C_flow = (np.pi / 4) * doi ** 2
        R_phpw_1h = doi / (2 * (doi + dio))
        R_phpw_2h = doi / (doi + dio)
        A_phpw_1h = dio / (2 * (doi + dio))
        A_phpw_2h = 1
        A_phpw_3h = doi / (doi + dio)

        if geo == "R":
            if hsur == 1:
                dtOSV = (
                    (qq * R_heated_1h)
                    / (R_flow * cpf * (1.25 * gsat * R_phpw_1h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            else:
                dtOSV = (
                    (qq * R_heated_2h)
                    / (R_flow * cpf * (1.25 * gsat * R_phpw_2h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        elif geo == "A":
            if hsur == 1:
                dtOSV = (
                    (qq * A_heated_1h)
                    / (A_flow * cpf * (1.25 * gsat * A_phpw_1h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            elif hsur == 2:
                dtOSV = (
                    (qq * A_heated_2h)
                    / (A_flow * cpf * (1.25 * gsat * A_phpw_2h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
            elif hsur == 2:
                dtOSV = (
                    (qq * A_heated_2h)
                    / (A_flow * cpf * (1.25 * gsat * A_phpw_3h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV,6), round(xOSV,6)
        elif geo == "C":
            dtOSV = (qq * C_heated) / (C_flow * cpf * (1.25 * gsat * (1.12 / p) ** 0.4)) / 100
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV,6), round(xOSV,6)

    def calOfiSZ(self, df):
        """
        :param df: dataframe
        :return: Heat flux [MW/m2]
        """
        # Time
        # 시작부분 코드

        start_time = time.time()
        # ------------------------
        # Saha and Zuber (1974) correlation
        print("Saha and Zuber Correlation (1974) OFI analysis")
        # Step 1
        # tsat     Saturation Temperature at P [K]
        # hfo      Saturated liquid Enthalpy at P [J/kg]
        # hgo      Saturated vapor Enthalpy at P [J/kg]
        # lam     Heat of vaporization at P [J/kg]
        # To Exit Temperature [K]

        for i in df.index:
            df.loc[i, "To"] = df.loc[i, "Ti"] + 1  # To 의 초기값설정

            for k in range(1, 10000):  # step 2 - 7 반복
                # Step 2
                # ho Enthalpy at exit [J/kg]
                df.loc[i, "ho"] = PropsSI(
                    "H", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )

                # Step 3
                # Pi Inlet pressure [bar]
                df.loc[i, "Pi_new"] = df.loc[i, "p"]

                for j in range(3):
                    df.loc[i, "Pi_old"] = df.loc[i, "Pi_new"]

                    # Step 4
                    # P_ave     Average pressure [bar]
                    # T_ave     Average Temperature [K]
                    # rho_ave   Density at Average pressure & Temperature [kg/m3]
                    # mu_ave    Viscosity at Average pressure & Temperature [Pa s]
                    # re_ave    reynolds number [-]

                    df.loc[i, "P_ave"] = 0.5 * (df.loc[i, "Pi_new"] + df.loc[i, "p"])
                    df.loc[i, "T_ave"] = 0.5 * (df.loc[i, "Ti"] + df.loc[i, "To"])
                    df.loc[i, "rho_ave"] = PropsSI(
                        "D",
                        "T",
                        df.loc[i, "T_ave"],
                        "P",
                        df.loc[i, "P_ave"] * 1e5,
                        df.loc[i, "refri"],
                    )
                    df.loc[i, "mu_ave"] = PropsSI(
                        "V",
                        "T",
                        df.loc[i, "T_ave"],
                        "P",
                        df.loc[i, "P_ave"] * 1e5,
                        df.loc[i, "refri"],
                    )
                    df.loc[i, "re_ave"] = PhysicalProperty.re(
                        self, df.loc[i, "g"], df.loc[i, "dh"], df.loc[i, "mu_ave"]
                    )

                    # f_ave Friction factor
                    if df.loc[i, "re_ave"] < 2300:
                        df.loc[i, "f_ave"] = 4 * (24 / df.loc[i, "re_ave"])
                    elif df.loc[i, "re_ave"] >= 4000:
                        df.loc[i, "f_ave"] = 4 * (
                            1.2810 ** -3 + 0.1143 * (df.loc[i, "re_ave"] ** (-0.3))
                        )
                    else:
                        df.loc[i, "f_ave"] = 4 * (
                            5.4 * 10 ** -3 + (2.3 * 1e-8) * (df.loc[i, "re_ave"] ** 0.75)
                        )

                    # Step 5
                    if df.loc[i, "flow"] == "Ho":  # in horizontal flow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (df.loc[i, "f_ave"] * df.loc[i, "lh"] / df.loc[i, "de"]) * 0.00001
                        )
                    elif df.loc[i, "flow"] == "Down":  # in vertical downflow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (
                                (df.loc[i, "f_ave"] * df.loc[i, "lh"] / df.loc[i, "de"])
                                - (df.loc[i, "rho_ave"] * 9.80665 * df.loc[i, "lh"])
                            )
                            * 0.00001
                        )
                    elif df.loc[i, "flow"] == "Up":  # in vertical upflow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (
                                df.loc[i, "f_ave"] * (df.loc[i, "lh"] / df.loc[i, "de"])
                                + (df.loc[i, "rho_ave"] * 9.80665 * df.loc[i, "lh"])
                            )
                            * 0.00001
                        )

                    df.loc[i, "P_diff"] = (
                        (df.loc[i, "Pi_old"] - df.loc[i, "Pi_new"]) * 100 / df.loc[i, "Pi_old"]
                    )

                    # Step 6
                # hi Enthalpy at inlet [J/kg]
                # q_guess Heat Flux based on heat balance [MW/m2]

                df.loc[i, "hi"] = PropsSI(
                    "H", "T", df.loc[i, "Ti"], "P", df.loc[i, "Pi_new"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "q_guess"] = PhysicalProperty.calQ(
                    self,
                    df.loc[i, "doi"],
                    df.loc[i, "dio"],
                    df.loc[i, "lh"],
                    df.loc[i, "g"],
                    df.loc[i, "geo"],
                    df.loc[i, "hsur"],
                    df.loc[i, "dh"],
                    df.loc[i, "ho"],
                    df.loc[i, "hi"],
                    df.loc[i, "de"],
                )

                # Step 7
                # Cpo       Cp at exit [J/kg K]
                # Ko        Thermal Conductivity at exit [W/m K]
                # Pe        Peclet number
                # Xo        Quality at exit

                df.loc[i, "Cpo"] = PropsSI(
                    "C", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "Ko"] = PropsSI(
                    "L", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "Pe"] = (
                    df.loc[i, "dh"] * df.loc[i, "g"] * df.loc[i, "Cpo"] / df.loc[i, "Ko"]
                )
                df.loc[i, "Xo"] = (df.loc[i, "ho"] - df.loc[i, "hfo"]) / (df.loc[i, "lam"])

                # Model
                # 1 SZ
                """
                Saha and Zuber (1974)
                """
                if df.loc[i, "Pe"] > 70000:
                    df.loc[i, "St_SZ"] = 0.0065
                else:
                    df.loc[i, "St_SZ"] = 455 / df.loc[i, "Pe"]

                # q_SZ      Heat Flux using the Saha-Zuber correlation [MW/m2]
                df.loc[i, "q_SZ"] = (
                    -df.loc[i, "g"]
                    * df.loc[i, "lam"]
                    * df.loc[i, "Xo"]
                    * df.loc[i, "St_SZ"]
                    * 10 ** -6
                )
                df.loc[i, "q_diff_SZ"] = df.loc[i, "q_SZ"] - df.loc[i, "q_guess"]

                # Error analysis
                if df.loc[i, "q_diff_SZ"] <= 0.01 and df.loc[i, "q_diff_SZ"] >= -0.01:
                    break
                elif df.loc[i, "q_diff_SZ"] > 0.01:
                    df.loc[i, "To"] = df.loc[i, "To"] + 0.1
                    continue
                elif df.loc[i, "q_diff_SZ"] < -0.01:
                    df.loc[i, "To"] = df.loc[i, "To"] - 0.1
                    continue
            print(
                "Saha and Zuber (1974): {:d}번째 데이터 완료, 최적화는 {:d}번 / Pi_new: {:.4f} / P_diff: {:.2f} / q_guess: {:.4f} / q_SZ: {:.4f} / q_diff_SZ: {:.2f}".format(
                    i + 1,
                    k,
                    df.loc[i, "Pi_new"],
                    df.loc[i, "P_diff"],
                    df.loc[i, "q_guess"],
                    df.loc[i, "q_SZ"],
                    df.loc[i, "q_diff_SZ"],
                )
            )

        # ----------------------------
        # 종료부분 코드
        print("start_time", start_time)
        print("--- %s seconds ---" % (time.time() - start_time))

    def calOfiJeong(self, df):
        # Jeong and Shim (2018) correlation
        print("Jeong and Shim Correlation (2019) OFI analysis")

        start_time = time.time()
        # ------------------------
        # Step 1
        # tsat     Saturation Temperature at P [K]
        # hfo      Saturated liquid Enthalpy at P [J/kg]
        # hgo      Saturated vapor Enthalpy at P [J/kg]
        # lam     Heat of vaporization at P [J/kg]
        # To Exit Temperature [K]

        for i in df.index:
            df.loc[i, "To"] = df.loc[i, "Ti"] + 1  # To 의 초기값설정

            for k in range(1, 10000):  # step 2 - 7 반복
                # Step 2
                # ho Enthalpy at exit [J/kg]
                df.loc[i, "ho"] = PropsSI(
                    "H", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )

                # Step 3
                # Pi Inlet pressure [bar]
                df.loc[i, "Pi_new"] = df.loc[i, "p"]

                for j in range(3):
                    df.loc[i, "Pi_old"] = df.loc[i, "Pi_new"]

                    # Step 4
                    # P_ave     Average pressure [bar]
                    # T_ave     Average Temperature [K]
                    # rho_ave   Density at Average pressure & Temperature [kg/m3]
                    # mu_ave    Viscosity at Average pressure & Temperature [Pa s]
                    # re_ave    reynolds number [-]

                    df.loc[i, "P_ave"] = 0.5 * (df.loc[i, "Pi_new"] + df.loc[i, "p"])
                    df.loc[i, "T_ave"] = 0.5 * (df.loc[i, "Ti"] + df.loc[i, "To"])
                    df.loc[i, "rho_ave"] = PropsSI(
                        "D",
                        "T",
                        df.loc[i, "T_ave"],
                        "P",
                        df.loc[i, "P_ave"] * 1e5,
                        df.loc[i, "refri"],
                    )
                    df.loc[i, "mu_ave"] = PropsSI(
                        "V",
                        "T",
                        df.loc[i, "T_ave"],
                        "P",
                        df.loc[i, "P_ave"] * 1e5,
                        df.loc[i, "refri"],
                    )
                    df.loc[i, "re_ave"] = PhysicalProperty.calRe(
                        self, df.loc[i, "g"], df.loc[i, "dh"], df.loc[i, "mu_ave"]
                    )

                    # f_ave Friction factor
                    if df.loc[i, "re_ave"] < 2300:
                        df.loc[i, "f_ave"] = 4 * (24 / df.loc[i, "re_ave"])
                    elif df.loc[i, "re_ave"] >= 4000:
                        df.loc[i, "f_ave"] = 4 * (
                            1.2810 ** -3 + 0.1143 * (df.loc[i, "re_ave"] ** (-0.3))
                        )
                    else:
                        df.loc[i, "f_ave"] = 4 * (
                            5.4 * 10 ** -3 + (2.3 * 1e-8) * (df.loc[i, "re_ave"] ** 0.75)
                        )

                    # Step 5
                    # Calculate Pi_new according to flow direction
                    if df.loc[i, "flow"] == "Ho":  # in horizontal flow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (df.loc[i, "f_ave"] * df.loc[i, "lh"] / df.loc[i, "de"]) * 0.00001
                        )
                    elif df.loc[i, "flow"] == "Down":  # in vertical downflow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (
                                (df.loc[i, "f_ave"] * df.loc[i, "lh"] / df.loc[i, "de"])
                                - (df.loc[i, "rho_ave"] * 9.80665 * df.loc[i, "lh"])
                            )
                            * 0.00001
                        )
                    elif df.loc[i, "flow"] == "Up":  # in vertical upflow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (
                                df.loc[i, "f_ave"] * (df.loc[i, "lh"] / df.loc[i, "de"])
                                + (df.loc[i, "rho_ave"] * 9.80665 * df.loc[i, "lh"])
                            )
                            * 0.00001
                        )

                    df.loc[i, "P_diff"] = (
                        (df.loc[i, "Pi_old"] - df.loc[i, "Pi_new"]) * 100 / df.loc[i, "Pi_old"]
                    )

                # Step 6
                # hi Enthalpy at inlet [J/kg]
                # q_guess Heat Flux based on heat balance [MW/m2]

                df.loc[i, "hi"] = PropsSI(
                    "H", "T", df.loc[i, "Ti"], "P", df.loc[i, "Pi_new"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "q_guess"] = PhysicalProperty.calQ(
                    self,
                    df.loc[i, "doi"],
                    df.loc[i, "dio"],
                    df.loc[i, "lh"],
                    df.loc[i, "g"],
                    df.loc[i, "geo"],
                    df.loc[i, "hsur"],
                    df.loc[i, "dh"],
                    df.loc[i, "ho"],
                    df.loc[i, "hi"],
                    df.loc[i, "de"],
                )

                # Step 7
                # Cpo       Cp at exit [J/kg K]
                # Ko        Thermal Conductivity at exit [W/m K]
                # Pe        Peclet number
                # Xo        Quality at exit

                df.loc[i, "Cpo"] = PropsSI(
                    "C", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "Ko"] = PropsSI(
                    "L", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "Pe"] = (
                    df.loc[i, "dh"] * df.loc[i, "g"] * df.loc[i, "Cpo"] / df.loc[i, "Ko"]
                )
                df.loc[i, "Xo"] = (df.loc[i, "ho"] - df.loc[i, "hfo"]) / (df.loc[i, "lam"])

                # Models
                # 1 Jeong and Shim
                """
                Jeong and Shim (2019)
                """
                if df.loc[i, "Pe"] > 70000:
                    df.loc[i, "St_JS"] = 0.15 / (df.loc[i, "Pe"] ** (0.25))
                else:
                    df.loc[i, "St_JS"] = 500 / df.loc[i, "Pe"]

                # q_JS      Heat Flux using the Jeong and Shim correlation [MW/m2]
                df.loc[i, "q_JS"] = (
                    -df.loc[i, "g"]
                    * df.loc[i, "lam"]
                    * df.loc[i, "Xo"]
                    * df.loc[i, "St_JS"]
                    * (10 ** -6)
                )
                df.loc[i, "q_diff_JS"] = df.loc[i, "q_JS"] - df.loc[i, "q_guess"]

                # Error analysis
                if df.loc[i, "q_diff_JS"] <= 0.01 and df.loc[i, "q_diff_JS"] >= -0.01:
                    break
                elif df.loc[i, "q_diff_JS"] > 0.01:
                    df.loc[i, "To"] = df.loc[i, "To"] + 0.1
                    continue
                elif df.loc[i, "q_diff_JS"] < -0.01:
                    df.loc[i, "To"] = df.loc[i, "To"] - 0.1
                    continue

            print(
                "Jeong and Shim (2019): {:d}번째 데이터 완료, 최적화는 {:.2f}번 / Pi_new: {:.4f} / P_diff: {:.2f} / q_guess: {:.4f} / q_JS: {:.4f} / q_diff_JS: {:.2f}".format(
                    i + 1,
                    k,
                    df.loc[i, "Pi_new"],
                    df.loc[i, "P_diff"],
                    df.loc[i, "q_guess"],
                    df.loc[i, "q_JS"],
                    df.loc[i, "q_diff_JS"],
                )
            )

        # ----------------------------
        # 종료부분 코드
        print("start_time", start_time)
        print("--- %s seconds ---" % (time.time() - start_time))

    def calOfiWF(self, df):
        # Whittle and Forgan (1967) correlation
        print("Whittle and Forgan Correlation (1967) OFI analysis")
        start_time = time.time()
        # ------------------------
        # Step 1
        # tsat     Saturation Temperature at P [K]
        # hfo      Saturated liquid Enthalpy at P [J/kg]
        # hgo      Saturated vapor Enthalpy at P [J/kg]
        # lam     Heat of vaporization at P [J/kg]
        # To Exit Temperature [K]

        for i in df.index:
            df.loc[i, "eta"] = 25
            # eta=14+(0.1*P(k)*0.986923)
            df.loc[i, "R"] = 1 / (1 + (df.loc[i, "eta"] * df.loc[i, "dh"] / df.loc[i, "lh"]))
            df.loc[i, "To"] = df.loc[i, "Ti"] + df.loc[i, "R"] * (
                df.loc[i, "tsat"] - df.loc[i, "Ti"]
            )

            for k in range(1, 5000):  # step 2 - 7 반복
                # Step 2
                # ho Enthalpy at exit [J/kg]
                df.loc[i, "ho"] = PropsSI(
                    "H", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )

                # Step 3
                # Pi Inlet pressure [bar]
                df.loc[i, "Pi_new"] = df.loc[i, "p"]

                for j in range(3):
                    df.loc[i, "Pi_old"] = df.loc[i, "Pi_new"]

                    # Step 4
                    # P_ave     Average pressure [bar]
                    # T_ave     Average Temperature [K]
                    # rho_ave   Density at Average pressure & Temperature [kg/m3]
                    # mu_ave    Viscosity at Average pressure & Temperature [Pa s]
                    # re_ave    reynolds number [-]

                    df.loc[i, "P_ave"] = 0.5 * (df.loc[i, "Pi_new"] + df.loc[i, "p"])
                    df.loc[i, "T_ave"] = 0.5 * (df.loc[i, "Ti"] + df.loc[i, "To"])
                    df.loc[i, "rho_ave"] = PropsSI(
                        "D",
                        "T",
                        df.loc[i, "T_ave"],
                        "P",
                        df.loc[i, "P_ave"] * 1e5,
                        df.loc[i, "refri"],
                    )
                    df.loc[i, "mu_ave"] = PropsSI(
                        "V",
                        "T",
                        df.loc[i, "T_ave"],
                        "P",
                        df.loc[i, "P_ave"] * 1e5,
                        df.loc[i, "refri"],
                    )
                    df.loc[i, "re_ave"] = PhysicalProperty.calRe(
                        self, df.loc[i, "g"], df.loc[i, "dh"], df.loc[i, "mu_ave"]
                    )

                    # f_ave Friction factor
                    if df.loc[i, "re_ave"] < 2300:
                        df.loc[i, "f_ave"] = 4 * (24 / df.loc[i, "re_ave"])
                    elif df.loc[i, "re_ave"] >= 4000:
                        df.loc[i, "f_ave"] = 4 * (
                            1.2810 ** -3 + 0.1143 * (df.loc[i, "re_ave"] ** (-0.3))
                        )
                    else:
                        df.loc[i, "f_ave"] = 4 * (
                            5.4 * 10 ** -3 + (2.3 * 1e-8) * (df.loc[i, "re_ave"] ** 0.75)
                        )

                    # Step 5
                    if df.loc[i, "flow"] == "Ho":  # in horizontal flow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (df.loc[i, "f_ave"] * df.loc[i, "lh"] / df.loc[i, "de"]) * 0.00001
                        )
                    elif df.loc[i, "flow"] == "Down":  # in vertical downflow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (
                                (df.loc[i, "f_ave"] * df.loc[i, "lh"] / df.loc[i, "de"])
                                - (df.loc[i, "rho_ave"] * 9.80665 * df.loc[i, "lh"])
                            )
                            * 0.00001
                        )
                    elif df.loc[i, "flow"] == "Up":  # in vertical upflow
                        df.loc[i, "Pi_new"] = (
                            df.loc[i, "Pi_old"]
                            + (
                                df.loc[i, "f_ave"] * (df.loc[i, "lh"] / df.loc[i, "de"])
                                + (df.loc[i, "rho_ave"] * 9.80665 * df.loc[i, "lh"])
                            )
                            * 0.00001
                        )

                    df.loc[i, "P_diff"] = (
                        (df.loc[i, "Pi_old"] - df.loc[i, "Pi_new"]) * 100 / df.loc[i, "Pi_old"]
                    )

                # Step 6
                # hi Enthalpy at inlet [J/kg]
                # q_guess Heat Flux based on heat balance [MW/m2]

                df.loc[i, "hi"] = PropsSI(
                    "H", "T", df.loc[i, "Ti"], "P", df.loc[i, "Pi_new"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "q_guess"] = PhysicalProperty.calQ(
                    self,
                    df.loc[i, "doi"],
                    df.loc[i, "dio"],
                    df.loc[i, "lh"],
                    df.loc[i, "g"],
                    df.loc[i, "geo"],
                    df.loc[i, "hsur"],
                    df.loc[i, "dh"],
                    df.loc[i, "ho"],
                    df.loc[i, "hi"],
                    df.loc[i, "de"],
                )

                # Step 7
                # Cpo       Cp at exit [J/kg K]
                # Ko        Thermal Conductivity at exit [W/m K]
                # Pe        Peclet number
                # Xo        Quality at exit

                df.loc[i, "Cpo"] = PropsSI(
                    "C", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "Ko"] = PropsSI(
                    "L", "T", df.loc[i, "To"], "P", df.loc[i, "p"] * 1e5, df.loc[i, "refri"]
                )
                df.loc[i, "Pe"] = (
                    df.loc[i, "dh"] * df.loc[i, "g"] * df.loc[i, "Cpo"] / df.loc[i, "Ko"]
                )
                df.loc[i, "Xo"] = (df.loc[i, "ho"] - df.loc[i, "hfo"]) / (df.loc[i, "lam"])

                # Models
                # 1 Whittle and Forgan
                """
                Whittle and Forgan (1967)
                """
                df.loc[i, "q_WF"] = (
                    (df.loc[i, "dh"] * df.loc[i, "g"] * (df.loc[i, "ho"] - df.loc[i, "hi"]))
                    / 4
                    / df.loc[i, "lh"]
                    * 10 ** -6
                )
                df.loc[i, "q_diff_WF"] = (df.loc[i, "q"] - df.loc[i, "q_WF"]) / df.loc[i, "q"]

                # Error analysis
                if df.loc[i, "q_diff_WF"] <= 0.1 and df.loc[i, "q_diff_JS"] >= -0.1:
                    break
                elif df.loc[i, "q_diff_WF"] > 0.1:
                    df.loc[i, "To"] = df.loc[i, "To"] + 0.1
                    continue
                elif df.loc[i, "q_diff_WF"] < -0.1:
                    df.loc[i, "To"] = df.loc[i, "To"] - 0.1
                    continue

            print(
                "Whittle and Forgan (1967): {:d}번째 데이터 완료, 최적화는 {:.2f}번 / Pi_new: {:.4f} / P_diff: {:.2f} / q_guess: {:.4f} / q_WF: {:.4f} / q_diff_WF: {:.2f}".format(
                    i + 1,
                    k,
                    df.loc[i, "Pi_new"],
                    df.loc[i, "P_diff"],
                    df.loc[i, "q_guess"],
                    df.loc[i, "q_WF"],
                    df.loc[i, "q_diff_WF"],
                )
            )

        # ----------------------------
        # 종료부분 코드
        print("start_time", start_time)
        print("--- %s seconds ---" % (time.time() - start_time))

    # statistic margin method
    def calMPE(self, y_true, y_pred):
        return np.mean((y_true - y_pred) / y_true)

    def calMSE(self, y_true, y_pred):
        return np.mean(np.square((y_true - y_pred)))

    def calMAPE(self, y_true, y_pred):
        return np.mean(np.abs((y_true - y_pred) / y_true)) * 100

    # The equations for CHF algorithm (Local hypothesis)
    def calCHFPark(self, rdcp, dh, g, xt_cal):
        """
        Park (2004) CHF correlation
        """
        # CHF 계산 (Park)
        alpha = round(0.71+4.6*rdcp-5.33*rdcp**2,16)
        gamma = round(0.1-0.58*rdcp+1.985*rdcp**2-1.54*rdcp**3,16)
        k1 = round(-0.343+0.22626*np.log(g)-0.01409*np.log(g)**2,16)
        k2 = 0.545
        k3 = round(2.6404-6.5*k1+6.1565*k1**2,16)
        fxt = round(np.sqrt(xt_cal*((1+xt_cal**2)**3)),16)
        q_cal = round((alpha/(dh ** k1)) * np.exp(-gamma*((g**k2)*fxt)**k3),16) # Park
        return alpha, gamma, k1, k2, k3, fxt, q_cal

    def calCHFDeng(self, rdcp=0.1, dh=0.1, g=0.1, xt_cal=0.5):
        """
        Deng (1997) CHF correlation
        """
        # CHF 계산 (Deng)
        alpha = round(1.669-6.544*(rdcp-0.448)**2,12)
        gamma = round(0.06523 + (0.1045/(math.sqrt((2*np.pi)*(math.log(rdcp))**2))) * math.exp(-5.413*((math.log(rdcp)+0.4537)**2/(math.log(rdcp)**2))),16)
        zxt = round((1+xt_cal**2)**3,12)
        #zxt = xt_cal*(1+math.exp(xt_cal))/math.exp(-xt_cal)
        q_cal = round((alpha/math.sqrt(dh)) * math.exp(-gamma*(g*xt_cal*zxt)**0.5),12) # Deng
        return alpha, gamma, zxt, q_cal

    def calCHFJeong(self, dh, g, cpv, cpf, rhov, rhof, lam, xt_cal=0.5):
        """
        Jeong (2023) CHF correlation
        """
        # calculate scaling factor
        S = (rhov*cpv/(cpf*rhof))**0.5
        # Set alpha parameters
        a_y0 = -0.08253
        a_xc = 0.3534
        a_A = 2.0986
        a_w1 = 0.4982
        a_w2 = 0.06853
        a_w3 = 0.15166

        # Set gamma parameters
        b_y0 = 0.0652
        b_xc = 0.472
        b_A = 0.1472
        b_w1 = 0.1172
        b_w2 = 0.0573
        b_w3 = 0.0266   
        zxt = (1+xt_cal**2)**3
        alpha = a_y0 + a_A*(1/(1+math.exp(-(S-a_xc+a_w1/2)/a_w2)))*(1-1/(1+math.exp(-(S-a_xc-a_w1/2)/a_w3)))
        gamma = b_y0 + b_A*(1/(1+math.exp(-(S-b_xc+b_w1/2)/b_w2)))*(1-1/(1+math.exp(-(S-b_xc-b_w1/2)/b_w3)))

        try:
            q_cal = round((alpha/math.sqrt(dh)) * math.exp(-gamma*(g*xt_cal*zxt)**0.5),12) # Deng
            return alpha, gamma, zxt, q_cal
        except:
            print("this step occurs an error: q_cal back to ")
            q_cal = round((alpha/math.sqrt(dh)) * math.exp(-gamma*(g*xt_cal*zxt)**0.5),12) # Deng
            return alpha, gamma, zxt, q_cal

    def sub_find_critical(self, dh, lh, g, q, lam, rdcp, Xi, Xe_ass, st_cal, modCHF, stepsize, tolerance):
        boolean=1
        # q_cal 계산
        if modCHF == 'Park':
            q_cal = lambda Xt: self.calCHFPark(rdcp, dh, g, Xt)[6]
        elif modCHF == 'Deng':
            q_cal = lambda Xt: self.calCHFDeng(rdcp, dh, g, Xt)[3]
        else:
            q_cal = q
        Xosv=lambda Xt: -q_cal(Xt)*10**6 / st_cal / g / lam
        Xe=lambda Xt: Xi + 4*q_cal(Xt)*10**6*lh / (lam*g*dh)  
        Xb=lambda Xt: max((Xosv(Xt),Xi))
        f=lambda Xt: Xosv(Xt)*np.log((Xe(Xt)-Xt)/Xb(Xt))+np.log((1-Xe(Xt)+Xosv(Xt)-Xosv(Xt)*Xt)/(1-Xb(Xt)+Xosv(Xt)))
        
        Critical_side=0
        Another_side=0

        if Xi<0:
            Xt = 0.1
            if Xe_ass < 0:
                old_Xt = 1e-6
            else:
                old_Xt = Xe_ass
            #Critical_side = 1e-6
            #Another_side = 0.5
            #converged = 'Converged'
            #return Critical_side, Another_side, converged
            
            Q=(Xosv(old_Xt)-Xosv(Xt))/(Xe(old_Xt)-Xe(Xt))
            Xe_over=(Q*Xe(Xt)-Xosv(Xt))/(Q-1)
            cnt_crt = 0

            while 1:
                # 함수 갱신텀
                fxt = f(Xt)
                xosv_tmp = Xosv(Xt)
                xe_tmp = Xe(Xt)

                #계산
                right_side=max([(1-xe_tmp+xosv_tmp)/xosv_tmp, xe_tmp, 0])   # 왼쪽 경계값을 결정  
                delta=abs(Xt-old_Xt)
                
                if right_side > Xt:
                    if boolean==1:
                        old_Xt=Xt
                        Xt=Xt+delta
                    else:
                        old_Xt=Xt
                        Xt=Xt+delta/2
                else:
                    if boolean==1 and (fxt>f(Xt+stepsize)) and fxt > 0:
                        old_Xt=Xt
                        Xt=Xt+delta/2
                    elif boolean==1:
                        boolean=0
                    else:
                        pass

                    if (fxt<0) and Another_side==0:
                        Another_side=Xt
                    else:
                        pass

                    if (fxt>f(Xt+stepsize) and fxt>0 and Critical_side==0) or (abs((old_Xt-Xt)/Xt) <= tolerance and Critical_side==0):
                        Critical_side=Xt
                        if Another_side==0:
                            old_Xt=Xt
                            Xt=Xt+delta
                            continue
                        else:
                            pass
                    else:
                        pass
                    old_Xt=Xt
                    Xt=Xt-delta/2
                cnt_crt +=1

                if cnt_crt > 5000:
                    converged = 'Diverged'
                    #print("Crtical_side : {}, Another_side : {}".format(Critical_side, Another_side))
                    return Critical_side, Another_side, converged
                else:
                    pass

                if (Critical_side != 0) and (Another_side != 0):
                    converged = 'Converged'
                    #print("Crtical_side : {}, Another_side : {}".format(Critical_side, Another_side))
                    return Critical_side, Another_side, converged
                else:
                    pass
        else:
            Xt= 1
            old_Xt=Xe_ass
            cnt_crt = 0
            
            while 1:
                # 함수 갱신 텀
                fxt = f(Xt)
                xosv_tmp = Xosv(Xt)
                xe_tmp = Xe(Xt)
                # 계산
                critical_min=max([(1-xe_tmp+xosv_tmp)/xosv_tmp, 0])    # 왼쪽 경계값을 결정                
                critical_max=xe_tmp
                delta=abs(Xt-old_Xt)
                if Xt > critical_max:
                    if boolean==1:
                        Xt=Xt-delta
                        old_Xt=Xt+delta/2
                    else:
                        old_Xt=Xt
                        Xt=Xt-delta/2
                else:
                    if boolean==1:
                        boolean=0
                    else:
                        pass

                    if (fxt>0 and Critical_side==0) or (abs((old_Xt-Xt)/Xt) <= tolerance and Critical_side==0):
                        Critical_side=Xt
                        if Another_side==0:
                            Xt=Xt-delta/2
                            cnt_crt += 1
                            continue
                        else:
                            pass
                    else:
                        pass

                    if Xt>critical_min and fxt<0 and Another_side==0:
                        Another_side=Xt
                    else:
                        pass

                    old_Xt=Xt
                    Xt=Xt+delta/2
                    cnt_crt += 1

                if cnt_crt > 5e3:
                    converged = 'Diverged'
                    return Critical_side, Another_side, converged
                else:
                    pass

                if (Critical_side != 0) and (Another_side != 0):
                    converged = 'Converged'
                    return Critical_side, Another_side, converged
                else:
                    continue

    def sub_bi(self, Xi, dh, lh, g, q, rdcp, lam, Xe_ass, st_cal, modCHF, stepsize, tolerance): 
        # q_cal 계산
        if modCHF == 'Park':
            q_cal = lambda Xt: self.calCHFPark(rdcp, dh, g, Xt)[6]
        elif modCHF == 'Deng':
            q_cal = lambda Xt: self.calCHFDeng(rdcp, dh, g, Xt)[3]
        else:
            q_cal = q
        Xosv=lambda Xt: -q_cal(Xt) / st_cal / g / lam * 10**6
        Xe=lambda Xt: Xi + 4*q_cal(Xt)*lh*10**6 / (lam*g*dh)
        Xb=lambda Xt: max((Xosv(Xt),Xi))
        f=lambda Xt: Xosv(Xt)*np.log((Xe(Xt)-Xt)/Xb(Xt))+np.log((1-Xe(Xt)+Xosv(Xt)-Xosv(Xt)*Xt)/(1-Xb(Xt)+Xosv(Xt)))

        if Xi==0:
            Xt_ass, converged=self.Xb_0( Xi, q, lh, lam, dh, g, rdcp, modCHF)
        else:
            Critical_side, Another_side, converged=self.sub_find_critical(dh, lh, g, q, lam, rdcp, Xi, Xe_ass, st_cal, modCHF, stepsize, tolerance)
            Xt_ass, converged=self.bisection(Critical_side, Another_side, Xi, g, q, rdcp, st_cal, lam, dh, lh, modCHF, stepsize, tolerance)
            
        # 최종 결과를 CHF 모델로 예측
        Xt_pre=Xt_ass
        print("Xt_pre {}".format(Xt_pre))
        # q_cal 계산
        if modCHF == 'Park':
            q_cal_all = lambda Xt: self.calCHFPark(rdcp, dh, g, Xt)
            alpha, gamma, k1, k2, k3, fxt, q_pre = q_cal_all(Xt_pre)
            Xosv=-q_pre / st_cal / g / lam * 10**6
            Xe_pre=Xi + 4*q_pre*lh*10**6 / (lam*g*dh)
            return round(Xosv,6), round(Xe_pre,6), round(q_pre,6), round(Xt_pre,6), round(alpha,6), round(gamma,6), round(k1,6), round(k2,6), round(k3,6), round(fxt,6), converged
        elif modCHF == 'Deng':
            q_cal_all = lambda Xt: self.calCHFDeng(rdcp, dh, g, Xt)
            alpha, gamma, fxt, q_pre = q_cal_all(Xt_pre)
            Xosv = - q_pre / st_cal / g / lam * 10**6
            Xe_pre=Xi + 4*q_pre*lh*10**6 / (lam*g*dh)
            return round(Xosv,6), round(Xe_pre,6), round(q_pre,6), round(Xt_pre,6), round(alpha,6), round(gamma,6), round(fxt,6), converged
        else:
            q_cal = q
        
    def Xb_0(self, Xi, q, lh, lam, dh, g, rdcp, modCHF):
        # q_cal 계산
        if modCHF == 'Park':
            q_cal = lambda Xt: self.calCHFPark(rdcp, dh, g, Xt)[6]
        elif modCHF == 'Deng':
            q_cal = lambda Xt: self.calCHFDeng(rdcp, dh, g, Xt)[3]
        else:
            q_cal = q
        RIGHT=1
        LEFT=1e-6
        Xe=lambda Xt: Xi + 4*q_cal(Xt)*lh*10**6 / (lam*g*dh)
        ged=lambda Xt: Xe(Xt)-Xt
        
        while np.abs(RIGHT-LEFT) > 0.001: # 이분법의 찾는범위가 eps값보다 작아지면 while문을 빠져나옴        
            CENTER=(LEFT+RIGHT)/2 # 현재 왼쪽값과 오른쪽값의 중간지점을 계산             
            g_r=ged(RIGHT) # 구간 오른쪽 경계값을 계산                    
            g_c=ged(CENTER) # 구간 왼쪽 경계값을 계산
            if g_r * g_c > 0: # 양 함수값의 부호판별
                RIGHT=CENTER # 같다면 다른구간 탐색
            else:
                LEFT=CENTER # 다르다면 현재중간값을 오른쪽값으로 놓는다.
        root=CENTER
        converged = 'Converged'
        return root, converged
        
    def bisection(self, Critical_Side, Another_Side, Xi, g, q, rdcp, st_cal, lam, dh, lh, modCHF, stepsize, tolerance):
        #Xb가 음수냐 양수냐에 따라 그래프의 모양이 달라진다.
        #Xb가 음수이면 해 Xt는 크리티컬 경계(왼쪽경계)의 오른쪽에 존재하고,
        #Xb가 양수이면 해 Xt는 크리티컬 경계(오른쪽경계)의 왼쪽에 존재한다.
        #해 Xt가 두 개가 존재할 경우도 있는데 이 경우 크리티컬에 가까운 값이 정해이다.
        #그러므로 이분법을 수행할 때 크리티컬 경계점~센터점에서 부호 판정을 수행해야 한다.
        # q_cal 계산
        if modCHF == 'Park':
            q_cal = lambda Xt: self.calCHFPark(rdcp, dh, g, Xt)[6]
        elif modCHF == 'Deng':
            q_cal = lambda Xt: self.calCHFDeng(rdcp, dh, g, Xt)[3]
        else:
            q_cal = q
        Xosv=lambda Xt: -q_cal(Xt) / st_cal / g / lam * 10**3
        Xe=lambda Xt: Xi + 4*q_cal(Xt)*lh*10**6 / (lam*g*dh)
        Xb=lambda Xt: max(Xosv(Xt),Xi)
        f=lambda Xt: Xosv(Xt)*np.log((Xe(Xt)-Xt)/Xb(Xt))+np.log((1-Xe(Xt)+Xosv(Xt)-Xosv(Xt)*Xt)/(1-Xb(Xt)+Xosv(Xt)))

        if f(Critical_Side)<0 and f(Another_Side)<0:
            root=Critical_Side
            converged = 'Converged'
            return root, converged
        else:
            pass

        cnt = 0
        while abs(Critical_Side-Another_Side) > tolerance: # 이분법의 찾는범위가 eps값보다 작아지면 while문을 빠져나옴        
            CENTER=(Critical_Side+Another_Side)/2   # 현재 왼쪽값과 오른쪽값의 중간지점을 계산         
            f_r=f(Critical_Side)   # 크리티컬에 가까운 경계쪽의 함수값을 계산
            f_c=f(CENTER)   # 구간 중간의 함수값을 계산
            cnt +=1
            if f_r * f_c > 0:    # 양 함수값의 부호판별
                Critical_Side=CENTER   # 같다면 다른구간 탐색
                #print("Critical_Side : {}".format(Critical_Side))
            else:
                Another_Side=CENTER   # 다르다면 현재 구간을 이등분하여 탐색
                #print("Another_Side : {}".format(Another_Side))
            
            cnt += 1
            if cnt > 5000:
                converged = 'Diverged'
                return Critical_Side, converged
            else:
                pass

        if cnt == 0:
            root = Critical_Side - Another_Side
            converged = 'Converged'
            return root, converged
        else:
            root=CENTER
            converged = 'Converged'
            return root, converged
    
    def calAlgCHFKim(self, i, rdcp, dh, lh, g, q, xi, xout, xt_cal_old, st_cal, lam, modCHF = 'Deng', stepsize = 0.0001, tolerance = 0.0001):
        # 여기에 김경동 선배 알고리즘 추가
        if modCHF == "Park":
            xosv_cal, xeq, q_cal, xt_cal_new, alpha, gamma, k1, k2, k3, Fxt, converged = self.sub_bi(xi, dh, lh, g, q, rdcp, lam, xout, st_cal, modCHF, stepsize, tolerance)
            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
            return round(xosv_cal,6), round(xeq,6), round(q_cal,6), round(xt_cal_new,6), round(alpha, 6), round(gamma,6), round(k1,6), round(k2,6), round(k3,6), round(Fxt, 6), converged
        else:
            xosv_cal, xeq, q_cal, xt_cal_new, alpha, gamma, Fxt, converged = self.sub_bi(xi, dh, lh, g, q, rdcp, lam, xout, st_cal, modCHF, stepsize, tolerance)
            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
            return round(xosv_cal,6), round(xeq,6), round(q_cal,6), round(xt_cal_new,6), round(alpha, 6), round(gamma,6), round(Fxt, 6), converged

    def calAlgCHFPark(self, i, rdcp, dh, lh, g, q, xi, xout, xt_cal_old, st_cal, lam, modCHF = 'Deng', stepsize = 0.0001, tolerance = 0.0001, init_flag = 0):
        # Flag 변수 설정
        cnt = 0 # Xt 계산 알고리즘 반복회수 계산
        cnt_nan = 0 # Xt_cal_new return 시도 실패 회수 계산
        
        while 1: # Xt_old와 Xt_new의 수렴여부 판단
            # xt_cal 초기값 (New 제약조건)
            if init_flag == -1:
                xt_cal_old = 1e-6
                init_flag = 1 # 재귀문 1회만 돌고 종료를 위해 값을 반환
            else:
                if xout < 0 or xi == 0:
                    xt_cal_old = 1e-6
                    init_flag = 1
                else:
                    xt_cal_old = xout
                    init_flag = -1 # Xe 초기값 설정이 오류가 났을 경우 처리 flag
            
            while 1: # Xt 찾기
                if cnt_nan == 5e3: # 종료 조건
                    if init_flag == -1:
                        return self.calAlgCHFPark(i, rdcp, dh, lh, g, q, xi, xout, xt_cal_old, st_cal, lam, modCHF = modCHF, stepsize = stepsize, tolerance = tolerance, init_flag = -1)
                    else:
                        pass
                    Fxt = 0
                    converged = 'Diverged'
                    print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                    # q_cal 계산
                    if modCHF == 'Park':
                        alpha, gamma, k1, k2, k3, zxt, q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                        return round(xosv_cal,6), round(xeq,6), round(q_cal,6), round(xt_cal,6), round(alpha,6), round(gamma,6), round(k1,6), round(k2,6), round(k3,6), round(Fxt,6), converged
                    elif modCHF == 'Deng':
                        alpha, gamma, zxt, q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        return round(xosv_cal,6), round(xeq,6), round(q_cal,6), round(xt_cal,6), round(alpha,6), round(gamma,6), round(Fxt,6), converged
                    else: # It is impossible process
                        q_cal = q 
                        return round(xt_cal_new,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)        
                else:
                    pass

                # q_cal 계산
                if modCHF == 'Park':
                    q_cal = self.calCHFPark(rdcp, dh, g, xt_cal_old)[6]
                elif modCHF == 'Deng':
                    q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal_old)[3]
                else:
                    q_cal = q
                
                # Xosv, Xeq, Xb 계산
                xeq = round(xi + (4*q_cal*10**6*lh)/(lam*g*dh),16)
                xosv_cal = round(-(q_cal/(st_cal*g*lam))*(10**6),16) # XOSV 계산값 (based on qCHF)
                xb = round(max(xi, xosv_cal), 16)
                                
                xosv_cal = round(-(q_cal/(st_cal*g*lam))*(10**6),16) # XOSV 계산값 (based on qCHF)
                if xi == 0: # xi가 0으로 들어올 경우, 
                    xb = round(xosv_cal,16)
                else:
                    xb = round(max(xi, xosv_cal), 16)
                 # Park 모델을 사용하기 위한 임시

                # Gauss-Seidel calculation 계산하기
                try:
                    cnt += 1
                    C1 = round(1 - xeq + xosv_cal - xosv_cal * xt_cal_old,16)
                    C2 = 1 - xb + xosv_cal

                    if C1/C2 > 0:
                        C3_tmp = round(C1/C2,16)
                    else:
                        C3_tmp = 1e-20
                    C3 = -np.log(C3_tmp)/xosv_cal
                    C4 = np.exp(C3)
                finally:
                    xt_cal = xeq - xb*C4
                    # New 제약조건
                    if xt_cal < 0: # xeq < 0이면 결국 q는 nan이 발생. xeq <0일 경우 강제로 1e-4로 결정
                        xt_cal = 1e-4
                    else:
                        pass

                    if xt_cal > 0:
                        if np.abs(xt_cal/xt_cal_old-1) < tolerance:
                            Fxt = 0
                            converged = 'Numerical, Converged'
                            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                            # q_cal 계산
                            if modCHF == 'Park':
                                alpha, gamma, k1, k2, k3, zxt, q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                                return round(xosv_cal,6), round(xeq,6), round(q_cal,6), round(xt_cal,6), round(alpha,6), round(gamma,6), round(k1,6), round(k2,6), round(k3,6), round(Fxt,6), converged
                            elif modCHF == 'Deng':
                                alpha, gamma, zxt, q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                                return round(xosv_cal,6), round(xeq,6), round(q_cal,6), round(xt_cal,6), round(alpha,6), round(gamma,6), round(Fxt,6), converged
                            else: # It is impossible process
                                q_cal = q 
                                return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                        else:
                            if np.abs(xt_cal/xt_cal_old-1) > 1:
                                xt_cal_new = 0.9999 * xt_cal_old  + 0.0001 * xt_cal
                                xt_cal_old = xt_cal_new
                                tmp = xt_cal
                                C4_old = C4
                                cnt_nan += 1
                                continue
                            else:
                                xt_cal_new = 0.99 * xt_cal_old  + 0.01 * xt_cal
                                xt_cal_old = xt_cal_new
                                tmp = xt_cal
                                C4_old = C4
                                cnt_nan += 1
                                continue                    
                    else:
                        xt_cal_old = 0.9999*xt_cal_old + 0.0001 * xt_cal
                        tmp = xt_cal
                        cnt_nan +=1
                        C4_old  = C4
                        continue

    def calIntgrXt(self, i, rdcp, dh, lh, g, q, xi, xout, xt_cal_old, st_cal, lam, modCHF = 0, stepsize = 0.0001, tolerance = 0.0001, flag_q = 1): # Measured CHF 기반 알고리즘

        # Lambda 함수 설정
        func = lambda x: xosv_cal * np.log((xeq -x)/xb) + np.log((1-xeq+xosv_cal-xosv_cal*x)/(1-xb+xosv_cal))
        funky = lambda x: xosv_cal*np.log((xeq-x)/xb)
        line = lambda x: -np.log((1-xeq+xosv_cal-xosv_cal*x)/(1-xb+xosv_cal))
        if flag_q == 1: # Measured qCHF를 기반으로 Xt 계산하는 알고리즘
            # 앞서 계산한 OSV 상관식으로 계산한 st이 삽입된 xosv_cal 계산
            xosv_cal = round(-(q * 10**6)/ (st_cal * g * lam),6) # XOSV 계산값 (based on qCHF)
            
            if xi == 0: # xi가 0으로 들어올 경우, 
                xb = round(xosv_cal,6)
            else:
                xb = round(max(xi, xosv_cal), 6)
            xeq = round(xi+ (4*q*10**6*lh)/(lam*g*dh),12)

            try:
                conv_tag = 0
                tmp = xosv_cal * np.log(xeq/xb) + np.log((1-xeq+xosv_cal)/(1-xb+xosv_cal))
            except ZeroDivisionError as e1:
                print("Error is occurs e1 : {}".format(e1))
                tmp = np.nan
            except ValueError as e2:
                print("Error is occurs e2 : {}".format(e2))
            except KeyError as e3:
                print("Error is occurs e3 : {}".format(e3))
            finally:
                # xt=0에서 Fxt가 존재하는가? 
                if np.isinf(tmp) == 1 or np.isnan(tmp) == 1:
                    xt_cal = round(xeq,6)
                else:
                    if xeq > 0:
                     xt_cal = round(xeq,6)
                    else:
                        funky = lambda x: xosv_cal*np.log((xeq-x)/xb)
                        line = lambda x: -np.log((1-xeq+xosv_cal-xosv_cal*x)/(1-xb+xosv_cal))
                        result = findIntersection(funky, line, 0)
                        xt_cal = result[0]
                        conv_tag = 1

            cnt = 0 # Xt 계산 알고리즘 반복회수 최대값
            cnt_nan = 0
            tmp_list = []

            while 1:
                try:
                    if conv_tag == 1: # 초기화 단계에서 수렴했으면 CHF 계산으로 넘어감
                        Fxt = 0
                        converged = 'Analytical, Converged'
                        xt_cal = xt_cal
                        # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q  
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6) 
                    else:
                        pass
                    cnt += 1 # 시작 계산횟수 1회
                    Fxt = xosv_cal * np.log((xeq -xt_cal)/xb) + np.log((1-xeq+xosv_cal-xosv_cal*xt_cal)/(1-xb+xosv_cal))
                except ZeroDivisionError as e1:
                    print("Error is occurs e1 : {}".format(e1))
                    #pass
                except ValueError as e2:
                    print("Error is occurs e2 : {}".format(e2))
                    #pass
                except KeyError as e3:
                    print("Error is occurs e3 : {}".format(e3))
                    #pass         
                finally:
                    pass

                    # 종료조건 1 : 반복횟수 100회 이상
                    if cnt > 5000:
                        converged = 'Diverged'
                        # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q  
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                    
                    # 종료조건 2: nan이 5번 이상이면 종료
                    if cnt_nan > 5:
                        print('해가 없으므로, 강제로 수렴시킵니다.')
                        Fxt = round(Fxt,6)
                        converged = 'Xt=Xeq, Converged'
                        # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q  
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                    
                    # 종료조건 3 : xi, Fxt 값에 따른 수렴 여부 파악
                    if xi > 0: # inlet quality > 0이면 그래프 개형은 왼쪽으로 발산
                        if np.abs(Fxt) < tolerance * 1e2:
                            Fxt = round(Fxt,6)
                            converged = 'Numerical, Converged'
                            xt_cal = xt_cal
                            # q_cal 계산
                            if modCHF == 'Park':
                                q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                            elif modCHF == 'Deng':
                                q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                            else:
                                q_cal = q  
                            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                            return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                        else:
                            if Fxt < 0:
                                cnt += 1
                                if np.isnan(Fxt) == 1:
                                    cnt_nan += 1
                                    xt_cal -= stepsize
                                    continue
                                else:                        
                                    # Bisection Method로 값 찾기
                                    xt_cal = round(bisect(func, xt_cal + stepsize, xt_cal)[0],6)
                                    Fxt = round(Fxt, 6)
                                    converged = 'Bisection, Converged'
                                    # q_cal 계산
                                    if modCHF == 'Park':
                                        q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                                    elif modCHF == 'Deng':
                                        q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                                    else:
                                        q_cal = q  
                                    print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                                    return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                            else:
                                cnt += 1
                                if np.isnan(Fxt) == 1:
                                    Fxt = 0
                                    cnt_nan += 1
                                    xt_cal -= stepsize
                                    continue
                                else:
                                    Fxt = round(Fxt,6)
                                    xt_cal -= stepsize
                                    continue
                    elif xi < 0: # inlet quality < 0이면 그래프 개형은 오른쪽으로 발산
                        if np.abs(Fxt) < tolerance * 1e2:
                            Fxt = round(Fxt, 6)
                            converged = 'Numerical, Converged'
                            # q_cal 계산
                            if modCHF == 'Park':
                                q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                            elif modCHF == 'Deng':
                                q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                            else:
                                q_cal = q  
                            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                            return round(xt_cal,6), converged, Fxt, round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                        else:
                            if Fxt < 0:
                                cnt += 1
                                if np.isnan(Fxt) == 1:
                                    Fxt = 0
                                    cnt_nan += 1
                                    xt_cal += stepsize
                                    continue
                                else:
                                    # Bisection Method로 값 찾기
                                    xt_cal = round(bisect(func, xt_cal-stepsize, xt_cal)[0],6)
                                    Fxt = round(Fxt, 6)
                                    converged = 'Bisection, Converged'
                                    # q_cal 계산
                                    if modCHF == 'Park':
                                        q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                                    elif modCHF == 'Deng':
                                        q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                                    else:
                                        q_cal = q                                    
                                    print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                                    return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                            else:
                                cnt += 1
                                if np.isnan(Fxt) == 1:
                                    Fxt = 0
                                    cnt_nan += 1
                                    xt_cal += stepsize
                                    continue
                                else:
                                    Fxt = round(Fxt,6)
                                    xt_cal += stepsize # Bisection Method
                                    continue
                    else:
                        Fxt = 0
                        xt_cal = xeq
                        converged = 'Xi==0, Converged'
                        # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)[6]
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                    
            cnt_nan = 0 # nan값을 삭제하고 초기화
        else:
            # Flag 변수 설정
            cnt = 0 # Xt 계산 알고리즘 반복회수 계산
            cnt_nan = 0 # Xt_cal_new return 시도 실패 회수 계산
            
            while 1: # Xt_old와 Xt_new의 수렴여부 판단
                #print("수렴 여부 판단 cnt_nan = {}, cnt = {}".format(cnt_nan, cnt))            
                # Xt_cal_old 계산
                #xt_cal_old = 0.0001

                while 1: # Xt 찾기
                    if cnt_nan == 1e5: # 종료 조건
                        xt_cal_new = xt_cal_new
                        Fxt = 0
                        converged = 'Diverged'
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal_new,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                    else:
                        pass

                    # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal_old)[6]
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal_old)
                        else:
                            q_cal = q
                    
                    # Xosv, Xeq, Xb 계산
                    xeq = round(xi+ (4*q_cal*10**6*lh)/(lam*g*dh),6)
                    xosv_cal = round(-(q_cal/(st_cal*g*lam))*(10**6),6) # XOSV 계산값 (based on qCHF)
                    xb = round(max(xi, xosv_cal), 6)
                    """
                    # New 제약조건
                    if xeq < 0: # xeq < 0이면 결국 q는 nan이 발생. xeq <0일 경우 강제로 1e-4로 결정
                        xeq = 1e-6
                    else:
                        pass
                    """
                    xosv_cal = round(-(q_cal/(st_cal*g*lam))*(10**6),6) # XOSV 계산값 (based on qCHF)
                    if xi == 0: # xi가 0으로 들어올 경우, 
                        xb = round(xosv_cal,6)
                    else:
                        xb = round(max(xi, xosv_cal), 6)
                     # Park 모델을 사용하기 위한 임시

                    # Gauss-Jeidel calculation 계산하기
                    try:
                        cnt += 1
                        C1 = round(1 - xeq + xosv_cal - xosv_cal * xt_cal_old,6)
                        C2 = 1 - xb + xosv_cal

                        if C1/C2 > 0:
                            C3_tmp = round(C1/C2,6)
                        else:
                            C3_tmp = 1e-20
                        C3 = -np.log(C3_tmp)/xosv_cal
                        C4 = np.exp(C3)
                                                
                        # Xt_cal_new 계산
                        # xt_cal_new = 0.99 * xt_cal_init + 0.01 * xt_cal

                        #print(q_cal, " + ", C4)
                    finally:                            
                        # Xt_cal_new 계산
                        xt_cal = xeq - xb*C4
                        xt_cal_new = 0.99*xt_cal_old + 0.01 * xt_cal
                        
                        if np.abs(xt_cal_new/xt_cal_old - 1) < tolerance:
                            xt_cal_new = xt_cal_new
                            Fxt = 0
                            converged = 'Numerical, Converged'
                            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                            return round(xt_cal_new,6), converged, round(Fxt,6), round(xosv_cal,6), round(xeq,6), round(q_cal,6)
                        else:                            
                            if xt_cal > 0:
                                xt_cal_old = xt_cal_new
                                C4_old = C4
                                cnt_nan += 1
                                continue
                            else:
                                xt_cal_old = 0.9999*xt_cal_old + 0.0001 * xt_cal
                                cnt_nan +=1
                                C4_old  = C4
                                continue
                            """
                            xt_cal_old = xt_cal_new # Park 모델을 사용하기 위한 임시
                            cnt_nan += 1
                            #print(i, " + ", cnt_nan)
                            continue                     
                        """
    
    def cal_xt(self, xi, xosv, org_xe, xe = 0.1):
        xb = round(max(xi, xosv), 12)

        # Define the equation
        x, y = sp.symbols('x y')

        eq = xosv* sp.log((xe-x)/xb) + sp.log((1-xe+xosv-xosv*x)/(1-xb+xosv))
        # New rate equation
        #Xt = 0 => must Xosv
        #eq = xosv* sp.log((xe-x)/xosv) + sp.log((1-xe+xosv-xosv*x)/(1-xosv+xosv))

        if xb >= 0.0:
            #print("Saturated flow boinling condition")
            try:
                sol = round(sp.nsolve(eq, (0, xe), solver='bisect'), 12)
            except:
                sol = round(xe, 12)
        else:
            #print("Subcooled flow boiling condition")
            if xe >= 1.0:
                sol = round(xe, 12)
            else:
                try:
                    sol = round(sp.nsolve(eq, (xe, 1), solver='bisect'),12)
                except:
                    sol = round(xe,12)
    
        return sol

    def cal_new(self, i, rdcp, dh, lh, g, q, xi, xe, rhof, cpf, kf, Pe, lam):
        """
        1) cal alpha, gamma
        2) cal qchf,i
        2) cal xosv
        3) cal xt_old
        4) cal qchf
        5) cal xt_new
        6) comp xt_old vs. xt_new
        7) Yes -> qchf
        8) No -> go to 4
        """
        tolerance = 0.1
        
       
        while 1:
             # Prepare: calculate alpha, gamma
            alpha = round(1.669-6.544*(rdcp-0.448)**2,16)
            gamma = round(0.06523 + (0.1045/(np.sqrt((2*np.pi)*(np.log(rdcp))**2))) * np.exp(-5.413*((np.log(rdcp)+0.4537)**2/(np.log(rdcp)**2))),16)

            # Step 1: calculate initial Xosv 
            old_xosv = self.cal_sz(q, rhof, dh, g, cpf, kf, Pe, lam)    

            # Step 2: calculate initial Xt
            old_xt = self.cal_xt(xi, old_xosv, xe)

            # Step 3: calculate initial qCHF
            old_qchf = self.calCHFDeng(rdcp, dh, g, old_xt)

            # Step 4: calculate new Xosv
            new_xosv = self.cal_sz(old_qchf, rhof, dh, g, cpf, kf, Pe, lam)

            # Step 5: calculate new Xt
            new_xt = self.cal_xt(xi, new_xosv, xe)

            # Step 5: calculate new qCHF
            new_qchf = self.calCHFDeng(rdcp, dh, g, new_xt)

            # Evaluation
            temp_eval = abs((new_xt - old_xt)/old_xt)
            if temp_eval <= tolerance:
                return new_qchf
            else:
                old_qchf = (old_qchf + new_qchf)/2

