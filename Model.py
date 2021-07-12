# Import Modules
import time
from CoolProp.CoolProp import PropsSI
from PhysicalProperty import *
from Numeric import *


class ModelOSV(PhysicalProperty):
    def __init__(self):
        print("Model_OSV is successfully started.")

    # Models or correlations
    def calGriffith(self, q, g, cpf, lam):  # Griffith et al. (1958) 1/10
        dtOSV = 5 * (q * 10 ** 6) / (g * cpf / 10)
        xOSV = -cpf * dtOSV / lam
        return round(dtOSV, 4), round(xOSV, 4)

    def calHancox(self, q, cpf, lam, kf, de, re, pr):  # Hancox and Nicoll (1971) 1/100배
        h = 0.4 * (re ** 0.662) * pr * (kf / de)
        dtOSV = (q * 10 ** 6) / h
        xOSV = -cpf * dtOSV / lam
        return round(dtOSV, 4), round(xOSV, 4)

    def calCosta(self, geo, q, v, cpf, lam):  # Costa (1967) # costa 는 10배 정도 크게 나옴
        if geo == "C":
            dtOSV = 1.8 * (q) / (np.sqrt(v / 100))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        if geo == "R":
            dtOSV = 1.28 * (q) / (np.sqrt(v / 100))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = np.nan
            xOSV = np.nan
            return round(dtOSV, 4), round(xOSV, 4)

    def calThom(self, q, g, lam, cpf, hfo):  # Thom (1966)
        dtOSV = 0.02 * hfo * (q * 10 ** 6) / (g * lam)
        xOSV = -cpf * dtOSV / lam
        return round(dtOSV, 4), round(xOSV, 4)

    """
    def calStaub(self): # Staub (1968)
        pass

    def calRogers(self,): # Rogers et al. (1987)
        theta = 30
        c1 = 2 + 3*np.cos(theta) - np.cos(theta)**3
        c2 = np.pi - theta + np.cos(theta)*np.sin(theta)
        c3 = np.sin(theta)*(np.cos(theta-10) - np.cos(theta+10))
        cs = 58 / (theta + 5) + 0.14
        rb = (3/(4*np.pi)) * ((c2 / c1) * cd * (u**2/9.8)) * (np.sqrt(1 + (8*np.pi**2/3)*(c1*c3/c2**2)*(cs/cd**2)*(9.8 * sigma / (rhof*u ** 4))) - 1)
        reb = rhof * u * (2*rb) / muf
        if reb < 20:
            cd = 24 / reb
        else:
            cd = 1.22
        cc = np.sqrt(1+)
        f = 0.046*re**(-0.2)
        tau = (0.046/8)*re**(-0.2)
        YB = (rhof/muf)*np.sqrt(tau/rhof) * (1+np.cos(theta)) * (3/(4*np.pi))*(c2/c1)*(cd*muf**2/9.8)*cc 
    
    def calJinghui(self,): # Jinghui and Rogers (1988)
        pass
    """

    def calHa2005(self, q, dh, kf, cpf, lam, pe):  # Ha et al. (2005)
        if pe < 52000:
            dtOSV = (1 / 918.5) * ((q * 10 ** 6) * dh / kf) * pe ** 0.08
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = 34.84 * ((q * 10 ** 6) * dh / kf) * (1 / pe) ** 0.876
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

    def calHa2018(self, rhof, rhov, lam, cpf, bo, v):
        ui = v / (1.18 * (9.8 * (rhof - rhov) / rhof ** 2) ** 0.25)
        if ui <= 1.55:
            dtOSV = 7.29 * (lam / cpf) * bo ** 0.8203
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = 32.94 * (lam / cpf) * bo ** 0.9016
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

    def calDix(self, kf, q, dh, cpf, lam, re, pr):  # Dix (1971)
        h = 0.036 * (kf / dh) * (re ** 0.8) * (pr ** 0.3)
        dtOSV = 0.00135 * (q * 10 ** 6 / h) * np.sqrt(re)
        xOSV = -cpf * dtOSV / lam
        return round(dtOSV, 4), round(xOSV, 4)

    def calKalitvianski(self, q, dh, kf, cpf, lam, pe):  # Kalitvianski (2000)
        if pe <= 36400:
            hosv = hsat + (5 / 455) * (((q * 10 ** 6) * dh * cpf) / kf)
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = (173.408 / 0.0065) * ((q * 10 ** 6) * dh / kf) * (1 / pe)
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

    def calSekoguchi(self, q, g, cpf, lam):  # Sekoguchi et al. (1980)
        dtOSV = 13.5 * (lam / cpf) * ((q * 10 ** 6) / (lam * g)) ** 0.65
        xOSV = -cpf * dtOSV / lam
        return round(dtOSV, 4), round(xOSV, 4)

    def calSahaZuber(self, q, rhof, dh, g, cpf, kf, Pe, lam):  # Saha and Zuber (1974)
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
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = 153.8 * (q * (10 ** 6)) / (g * cpf)
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

    def calParkSahaZuber(self, q, rhof, dh, g, cpf, kf, Pe, lam):  # Park (2004)
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
            return round(dtOSV, 4), round(xOSV, 4)
        elif Pe > 200000:
            dtOSV = (
                0.08923 * np.exp(-(Pe * (kf ** 0.45) / (dh ** 0.53 * g ** 0.37)) / 25313.63287)
                + 0.00659 * np.exp(-(Pe * (kf ** 0.45) / (dh ** 0.53 * g ** 0.37)) / 211422.70151)
                + 0.00146
            )
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = 153.8 * (q * (10 ** 6)) / (g * cpf)
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

    def calMSZ(
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
                    return round(dtOSV, 4), round(xOSV, 4)
                else:
                    dtOSV = ((153.8 * (qq)) / (g * cpf)) * R_phpw_2h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV, 4), round(xOSV, 4)
            elif geo == "A":
                if hsur == 1:
                    dtOSV = ((153.8 * (qq)) / (g * cpf)) * A_phpw_1h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV, 4), round(xOSV, 4)
                else:
                    dtOSV = ((153.8 * (qq)) / (g * cpf)) * A_phpw_2h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV, 4), round(xOSV, 4)
            else:
                dtOSV = 153.8 * (qq) / (g * cpf)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        else:
            if geo == "R":
                if hsur == 1:
                    dtOSV = (0.0022 * (qq * dh) / kf) * R_phpw_1h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV, 4), round(xOSV, 4)
                else:
                    dtOSV = (0.0022 * (qq * dh) / kf) * R_phpw_2h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV, 4), round(xOSV, 4)
            elif geo == "A":
                if hsur == 1:
                    dtOSV = (0.0022 * (qq ** dh) / kf) * A_phpw_1h
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV, 4), round(xOSV, 4)
                else:
                    dtOSV = 0.0022 * (qq * dh) / kf
                    xOSV = -cpf * dtOSV / lam
                    return round(dtOSV, 4), round(xOSV, 4)
            else:
                dtOSV = 0.0022 * (qq * dh) / kf
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)

    def calLevy(self, sigma, dh, rhof, muf, kf, re, pr, cpf, g, q, lam, v):  # Levy (1967)
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
            return round(dtOSV, 4), round(xOSV, 4)
        elif YB <= 5:
            dtOSV = (q * 10 ** 6) / h - Q * pr * YB
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = (q * 10 ** 6) / h - 5 * Q * (pr + np.log(1 + pr * ((YB / 5) - 1)))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

    def calBowring(self, p, q, v, lam, cpf):  # Bowring (1960)
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
        return round(dtOSV, 4), round(xOSV, 4)

    def calUnal(self, q, pr, dh, v, cpf, kf, re, refri, lam):  # Unal (1975)
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
                return round(dtOSV, 4), round(xOSV, 4)
            elif v <= 0.45:
                dtOSV = (0.11 * (q * 10 ** 6)) / ((kf / dh) * 0.023 * re ** 0.8 * pr ** 0.4)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        else:
            if v > 0.45:
                dtOSV = (0.28 * (q * 10 ** 6)) / ((kf / dh) * 0.023 * re ** 0.8 * pr ** 0.4)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
            elif v <= 0.45:
                dtOSV = (0.11 * (q * 10 ** 6)) / ((kf / dh) * 0.023 * re ** 0.8 * pr ** 0.4)
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)

    def calEl(self, Bo_el, pr, lh, de, dtin, cpf, lam):  # El-Morshedy (2012)
        """
        :param Bo_el: Boiling number of El-Mosherdy [-]
        :param pr: prandtl number [-]
        :param lh: Heated length of channel [m]
        :param dh: Hydraulic diameter [m][
        :param dtin: Inlet liquid subcooling [K]
        :param cpf: Liquid specific heat [J/kg-K]
        :param lam: Heat of vaporization [J/kg]
        :return: Equilibrium thermal quality [-]
        """
        if dtin is None:
            dtOSV = np.nan
            xOSV = np.nan
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = dtin * (Bo_el ** 0.0094) * (pr ** 1.606) / ((lh / de) ** 0.533)
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

    def calJeong(self, q, rhof, rhov, dh, v, cpf, kf, Pe, lam, Pr, We, Ca):  # Jeong and Shim (2021)
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
        if We <= 200:
            dtOSV = (q * (10 ** 6)) / (rhof * v * cpf * (17.25 * Pe ** -0.75 * Ca ** -0.15))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)
        else:
            # dtOSV = (q * 10 ** 6 * dh) / (kf * (0.125 * Pe ** 0.75)) # for Nu
            dtOSV = (q * (10 ** 6)) / (
                rhof * v * cpf * (0.0965 * Pe ** -0.225 * (rhof / (rhof - rhov) ) ** 1.5)
            )  # for St
            # dtOSV = (q * (10 ** 6)) / (rhof * v * cpf * (0.0065)) # for St
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

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
                return round(dtOSV, 4), round(xOSV, 4)
            else:
                dtOSV = (qq * R_heated_2h) / (R_flow * cpf * ((gsat + 27) / 0.58))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        elif geo == "A":
            if hsur == 1:
                dtOSV = (qq * A_heated_1h) / (A_flow * cpf * ((gsat + 27) / 0.58))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
            else:
                dtOSV = (qq * A_heated_2h) / (A_flow * cpf * ((gsat + 27) / 0.58))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = (qq * C_heated) / (C_flow * cpf * ((gsat + 27) / 0.58))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

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
                return round(dtOSV, 4), round(xOSV, 4)
            else:
                dtOSV = (qq * R_heated_2h) / (R_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        elif geo == "A":
            if hsur == 1:
                dtOSV = (qq * A_heated_1h) / (A_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
            elif hsur == 2:
                dtOSV = (qq * A_heated_2h) / (A_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
            elif hsur == 3:
                dtOSV = (qq * A_heated_2h) / (A_flow * cpf * (gsat * 1.11))
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        else:
            dtOSV = (qq * C_heated) / (C_flow * cpf * (gsat * 1.11))
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

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
                return round(dtOSV, 4), round(xOSV, 4)
            else:
                dtOSV = (
                    (qq * R_heated_2h)
                    / (R_flow * cpf * (1.25 * gsat * R_phpw_2h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        elif geo == "A":
            if hsur == 1:
                dtOSV = (
                    (qq * A_heated_1h)
                    / (A_flow * cpf * (1.25 * gsat * A_phpw_1h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
            elif hsur == 2:
                dtOSV = (
                    (qq * A_heated_2h)
                    / (A_flow * cpf * (1.25 * gsat * A_phpw_2h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
            elif hsur == 2:
                dtOSV = (
                    (qq * A_heated_2h)
                    / (A_flow * cpf * (1.25 * gsat * A_phpw_3h * (1.12 / p) ** 0.4))
                    / 100
                )
                xOSV = -cpf * dtOSV / lam
                return round(dtOSV, 4), round(xOSV, 4)
        elif geo == "C":
            dtOSV = (qq * C_heated) / (C_flow * cpf * (1.25 * gsat * (1.12 / p) ** 0.4)) / 100
            xOSV = -cpf * dtOSV / lam
            return round(dtOSV, 4), round(xOSV, 4)

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
        alpha_park = round(0.71+4.6*rdcp-5.33*rdcp**2,6)
        gamma_park = round(0.1-0.58*rdcp+1.9851*rdcp**2-1.54*rdcp**3,6)
        k1_park = round(-0.343+0.22626*np.log(g)-0.01409*np.log(g)**2,6)
        k2_park = 0.545
        k3_park = round(2.6404-6.5*k1_park+6.1565*k1_park**2,6)
        fxt_park = round(np.sqrt(xt_cal*(1+xt_cal**2)**3),6)
        q_cal_park = round(alpha_park/(dh ** k1_park) * np.exp(-gamma_park*((g**k2_park)*fxt_park)**k3_park),6) # Park
        return q_cal_park

    def calCHFDeng(self, rdcp, dh, g, xt_cal):
        """
        Deng () CHF correlation
        """
        # CHF 계산 (Park)
        alpha_deng = round(1.669-6.544*(rdcp-0.448)**2,6)
        gamma_deng = round(0.06523 + (0.1045/(np.sqrt((2*np.pi)*(np.log(rdcp))**2))) * np.exp(-5.413*((np.log(rdcp)+0.4537)**2/(np.log(rdcp)**2))),6)
        fxt_deng = round(np.sqrt(g*xt_cal*((1+xt_cal**2)**3)),6)
        q_cal_park = round((alpha_deng/np.sqrt(dh)) * np.exp(-gamma_deng*fxt_deng),6) # Park
        return q_cal_park

    

    def calIntgrXt(self, i, rdcp, dh, lh, g, q, xi, xout, xt_cal_old, st_cal, lam, modCHF = 'Deng', stepsize = 0.0001, tolerance = 0.0001, flag_q = 1):
        # Lambda 함수 설정
        func = lambda x: xosv_cal * np.log((xeq -x)/xb) + np.log((1-xeq+xosv_cal-xosv_cal*x)/(1-xb+xosv_cal))
        funky = lambda x: xosv_cal*np.log((xeq-x)/xb)
        line = lambda x: -np.log((1-xeq+xosv_cal-xosv_cal*x)/(1-xb+xosv_cal))

        if flag_q == 0: # Measured qCHF를 기반으로 Xt 계산하는 알고리즘
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
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q  
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4) 
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
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q  
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
                    
                    # 종료조건 2: nan이 5번 이상이면 종료
                    if cnt_nan > 5:
                        print('해가 없으므로, 강제로 수렴시킵니다.')
                        Fxt = round(Fxt,6)
                        converged = 'Xt=Xeq, Converged'
                        # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q  
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
                    
                    # 종료조건 3 : xi, Fxt 값에 따른 수렴 여부 파악
                    if xi > 0: # inlet quality > 0이면 그래프 개형은 왼쪽으로 발산
                        if np.abs(Fxt) < tolerance * 1e2:
                            Fxt = round(Fxt,6)
                            converged = 'Numerical, Converged'
                            xt_cal = xt_cal
                            # q_cal 계산
                            if modCHF == 'Park':
                                q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                            elif modCHF == 'Deng':
                                q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                            else:
                                q_cal = q  
                            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                            return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
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
                                        q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                                    elif modCHF == 'Deng':
                                        q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                                    else:
                                        q_cal = q  
                                    print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                                    return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
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
                                q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                            elif modCHF == 'Deng':
                                q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                            else:
                                q_cal = q  
                            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                            return round(xt_cal,6), converged, Fxt, round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
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
                                        q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                                    elif modCHF == 'Deng':
                                        q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                                    else:
                                        q_cal = q                                    
                                    print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                                    return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
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
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal)
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal)
                        else:
                            q_cal = q
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
                    
            cnt_nan = 0 # nan값을 삭제하고 초기화
        elif flag_q == 1:        
            # Flag 변수 설정
            cnt = 0 # Xt 계산 알고리즘 반복회수 계산
            cnt_nan = 0 # Xt_cal_new return 시도 실패 회수 계산
            
            while 1: # Xt_old와 Xt_new의 수렴여부 판단
                #print("수렴 여부 판단 cnt_nan = {}, cnt = {}".format(cnt_nan, cnt))
                if cnt_nan == 1000:
                    converged = 'Diverged'
                    print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                    return round(xt_cal_new,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
                else:
                    while 1: # Xt 찾기
                        # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal_old)
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal_old)
                        else:
                            q_cal = q

                        # Xosv, Xeq, Xb 계산
                        xosv_cal = round(-(q_cal * 10**6)/ (st_cal * g * lam),6) # XOSV 계산값 (based on qCHF)
                        
                        if xi == 0: # xi가 0으로 들어올 경우, 
                            xb = round(xosv_cal,6)
                        else:
                            xb = round(max(xi, xosv_cal), 6)

                        xeq = round(xi+ (4*q_cal*10**6*lh)/(lam*g*dh),6)

                        # Xt 계산하기
                        try:
                            cnt += 1
                            xt_cal_init = xt_cal_old
                            while 1:
                                try: 
                                    Fxt_old = xosv_cal * np.log((xeq - xt_cal_old)/xb) + np.log((1-xeq+xosv_cal-xosv_cal*xt_cal_old)/(1-xb+xosv_cal)) 
                                finally:
                                    if np.isinf(Fxt_old) == 1 or np.isnan(Fxt_old) == 1: # 만약 Fxt가 없으면 xt_cal_old를 stepsize만큼 이동
                                        if xt_cal_old < 0 or xeq > 1: # Xt_cal_old, xeq의 값이 없으면 xout = xt와 같다.
                                            xt_cal_old = xeq
                                            Fxt = 1
                                            Fxt_old = 1
                                            break
                                        
                                        if xi >0:
                                            xt_cal_old -= stepsize
                                            continue
                                        else:
                                            xt_cal_old += stepsize
                                            continue
                                    else:
                                        if Fxt_old < 0:
                                            Fxt = 0
                                            xt_cal_old = xeq
                                            tmp = -1
                                            break
                                        else:
                                            Fxt = Fxt_old
                                            tmp = 1
                                            break  
                        finally:                            
                            # 종료조건 1 : 반복횟수 100회 이상
                            if cnt > 5000:
                                xt_cal = xt_cal_old
                                Fxt = Fxt_old
                                converged = 'Diverged'
                                break

                            # 종료조건 2 : xi, Fxt 값에 따른 수렴 여부 파악
                            if xi > 0: # inlet quality > 0이면 그래프 개형은 왼쪽으로 발산
                                if np.abs(Fxt/Fxt_old-1) < tolerance*10:
                                    Fxt = Fxt_old
                                    xt_cal = xt_cal_old
                                    converged = 'Numerical, Converged'
                                    break
                                else:
                                    if Fxt_old < 0 :
                                        if tmp == -1:
                                            Fxt = 0
                                            xt_cal = xt_cal_old
                                            converged = 'Numerical, Converged'
                                            break
                                        else:
                                            # Bisection Method로 값 찾기
                                            xt_cal = round(bisect(func, xt_cal_old+stepsize, xt_cal_old)[0],6)
                                            Fxt = 0
                                            converged = 'Bisection, Converged'
                                            break
                                    else:
                                        cnt += 1 # 계산회수 1회 추가
                                        Fxt = Fxt_old
                                        xt_cal_old -= stepsize
                                        continue
                            elif xi < 0: # inlet quality < 0이면 그래프 개형은 오른쪽으로 발산
                                if np.abs(Fxt/Fxt_old-1) < tolerance*10:
                                    Fxt = Fxt_old
                                    xt_cal = xt_cal_old
                                    converged = 'Numerical, Converged'
                                    break
                                else:
                                    if Fxt_old < 0 :
                                        if tmp == -1:
                                            Fxt = 0
                                            xt_cal = xt_cal_old
                                            converged = 'Numerical, Converged'
                                            break
                                        else:
                                            # Bisection Method로 값 찾기
                                            xt_cal = round(bisect(func, xt_cal_old+stepsize, xt_cal_old)[0],6)
                                            Fxt = 0
                                            converged = 'Bisection, Converged'
                                            break
                                    else:
                                        cnt += 1 # 계산회수 1회 추가
                                        Fxt = Fxt_old
                                        xt_cal_old += stepsize
                                        continue
                            else:
                                Fxt = 0
                                xt_cal = xeq
                                converged = 'Xi==0, Converged'
                                break
                # Xt_old Xt_new 계산
                
                xt_cal_new = 0.99 * xt_cal_init + 0.01 * xt_cal
                
                if np.abs(xt_cal_new - xt_cal_init) < tolerance:
                    print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                    return round(xt_cal_new,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
                else:
                    xt_cal_old = xt_cal_new
                    cnt_nan += 1
                    del Fxt_old
                    continue
        else:
            # Flag 변수 설정
            cnt = 0 # Xt 계산 알고리즘 반복회수 계산
            cnt_nan = 0 # Xt_cal_new return 시도 실패 회수 계산
            
            while 1: # Xt_old와 Xt_new의 수렴여부 판단
                #print("수렴 여부 판단 cnt_nan = {}, cnt = {}".format(cnt_nan, cnt))            
                # Xt_cal_old 계산
                xt_cal_old = 0.0001

                while 1: # Xt 찾기
                    if cnt_nan == 1e4: # 종료 조건
                        xt_cal_new = xt_cal_new
                        Fxt = 0
                        converged = 'Diverged'
                        print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                        return round(xt_cal_new,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
                    else:
                        pass

                    # q_cal 계산
                        if modCHF == 'Park':
                            q_cal = self.calCHFPark(rdcp, dh, g, xt_cal_old)
                        elif modCHF == 'Deng':
                            q_cal = self.calCHFDeng(rdcp, dh, g, xt_cal_old)
                        else:
                            q_cal = q
                    
                    # Xosv, Xeq, Xb 계산
                    xeq = round(xi+ (4*q_cal*10**6*lh)/(lam*g*dh),6)

                    if xeq < 0: # xeq < 0이면 결국 q는 nan이 발생. xeq <0일 경우 강제로 1e-4로 결정
                        xeq = 1e-6
                    else:
                        pass

                    xosv_cal = round(-(q_cal/(st_cal*g*lam))*(10**6),6) # XOSV 계산값 (based on qCHF)
                    if xi == 0: # xi가 0으로 들어올 경우, 
                        xb = round(xosv_cal,6)
                    else:
                        xb = round(max(xi, xosv_cal), 6)
                    
                    # Gauss-Jeidel calculation 계산하기
                    try:
                        cnt += 1
                        C1 = 1 - xeq + xosv_cal - xosv_cal * xt_cal_old
                        C2 = 1 - xb + xosv_cal

                        if C1/C2 >0:
                            C3 = -np.log(C1/C2)
                        else:
                            C3 = -np.log(10)
                        
                        C4 = np.exp(C3)
                        """
                        if cnt_nan == 0:
                            C4_old = 1
                            if C4/C4_old < 0:
                                C4 = 0.99*C4 + 0.01*C4_old
                            else:
                                pass
                        else:
                            if C4/C4_old < 0:
                                C4 = 0.99*C4 + 0.01*C4_old
                            else:
                                pass
                        """
                        # Xt_cal_new 계산
                        xt_cal = xeq - xb*C4
                        xt_cal_new = 0.99*xt_cal_old + 0.01 * xt_cal
                        #print(q_cal, " + ", C4)
                    finally:                            
                        #if np.abs(xt_cal_new-xt_cal_old) < 0.001:
                        if np.abs(xt_cal/xt_cal_old - 1) < 0.001:
                            xt_cal_new = xt_cal_new
                            Fxt = 0
                            converged = 'Numerical, Converged'
                            print("{}번째 데이터는 {}되었습니다.".format(i+1, converged))
                            return round(xt_cal_new,6), converged, round(Fxt,6), round(xosv_cal, 4), round(xeq, 4), round(q_cal, 4)
                        else:

                            if xt_cal > 0:
                                xt_cal_old = xt_cal_new
                                C4_old = C4
                            else:
                                xt_cal_old = 0.9999*xt_cal_old + 0.0001 * xt_cal
                                C4_old  = C4
                                
                            cnt_nan += 1
                            #print(i, " + ", cnt_nan)
                            continue                          


    