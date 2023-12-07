import numpy as np
import pandas as pd

class PhysicalProperty():
    def __init__(self):
        print("Property_OSV is successfully started.")

    # Dimensionless number
    def cal_pe(self, dh, g, cpf, kf):  # Peclet number
        Pe = (g * dh * cpf) / kf
        return Pe

    def cal_st(self, q, cpf, rhof, v, tosv):  # Stanton number
        St = (q * (10 ** 6)) / (cpf * rhof * v * tosv)
        return St

    def cal_we(self, dh, v, rhof, sigma):  # Weber number
        We = rhof * (v ** 2) * dh / sigma
        return We

    def cal_re(self, dh, g, muf):  # Reynolds number
        Re = (g * dh) / muf
        return Re

    def cal_bd(self, rhof, rhov, dh, sigma):  # Bond number
        Bd = (9.8 * (rhof - rhov) * (dh ** 2)) / (sigma)
        return Bd

    def cal_gr(self, dh, rhof, tosv, muf):
        pass
        #Gr = 9.8*(dh**3)*(rhof**2)*()
        #return Gr

    def cal_bo(self,q,lam,g): # Boiling number
        Bo = q*10**6/(lam*g)
        return Bo

    def cal_ga(self,rhof, dh, muf): # Galileo number == archimedes number
        Ar = (9.8 * rhof * dh ** 3) / muf ** 2
        return Ar

    def cal_jh(self, q, cpf, rhof, v, tosv, muf, kf ): # Colburn j factor (jH)
        Jh = ((q * (10 ** 6)) / (cpf * rhof * v * tosv))* ((cpf * muf) / kf)**(2/3)
        return Jh

    def cal_gz(self, dh, lh, re, pr): # Graetz number
        if lh is None:
            Gz = np.nan
        else:
            Gz = (dh/lh) * re * pr
        return Gz

    def cal_ja(self, cpf, tosv, lam): # Jakob number
        Ja = cpf * tosv / lam
        return Ja

    def cal_Z(self, muf, rhof, dh, sigma): # Ohnesorge number
        Z = muf / (rhof * dh * sigma) ** 0.5
        return Z

    def cal_fr(self, v, dh):  # Lewis number or Froude number
        Fr = (v ** 2) / (9.8 * dh)
        return Fr

    def cal_ca(self, muf, v, sigma, rhof): # capillary number
        Ca = (muf*v)/(sigma)
        return Ca

    def cal_co(self, tosv, rhof, lam, dh, kf, muf): # Condensation number
        Co = (9.8*(rhof**2)*lam*dh**3)/(kf*muf*tosv)
        return Co

    def cal_pr(self, cpf, muf, kf):  # Prandtl number
        Pr = (cpf * muf) / kf
        return Pr

    def cal_nu(self, q, dh, kf, tosv):  # Nusselt number
        Nu = ((q * (10 ** 6)) * dh) / (kf * tosv)
        return Nu

    def cal_ti(self, tsat, dtin):
        tin = tsat - 273.15 - dtin
        return tin

    def cal_qrt(self, q, doi, dio, geo, hs, g, cpf, dtin, lh):
        qq = q * (10 ** 6)
        R_heated_1h = doi * lh
        R_heated_2h = 2 * doi * lh
        R_flow = doi * dio
        A_heated_1h = (np.pi * dio) * lh
        A_heated_2h = (np.pi * (dio + doi)) * lh
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        C_heated = (np.pi) * doi * lh
        C_flow = (np.pi / 4) * doi ** 2

        if geo == 'R':
            if hs == 1:
                Qratio = (qq * R_heated_1h) / (R_flow * g * cpf * dtin)
                return Qratio
            else:
                Qratio = (qq * R_heated_2h) / (R_flow * g * cpf * dtin)
                return Qratio
        elif geo == 'A':
            if hs == 1:
                Qratio = (qq * A_heated_1h) / (A_flow * g * cpf * dtin)
                return Qratio
            else:
                Qratio = (qq * A_heated_2h) / (A_flow * g * cpf * dtin)
                return Qratio
        else:
            Qratio = (qq * C_heated) / (C_flow * g * cpf * dtin)
            return Qratio
        
    def cal_xe(self, qav, q, q2, doi, dio, geo, hs, g, xi, lh, lam, ch = 1):
        # calculate heated section parameter
        if hs == 1:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * dio)
            R_heated = doi
            qsol = qav * (10 ** 3) # Lambda = kJ/kg
        elif hs == 2:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * doi)
            R_heated = 2 * doi
            qsol = (qav) * (10 ** 3)
        else:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * (dio + doi))
            R_heated = 2 * doi
            qsol = qav * (10 ** 3)
        
        # calculate flow area
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        R_flow = doi * dio
        C_flow = (np.pi / 4) * doi ** 2

        # calculate quality
        if geo == 'R':
            xe = xi + (qsol * R_heated * lh) / (R_flow * g * lam)
            return xe                
        elif geo == 'A':
            if hs == 1:
                xe = xi + (np.pi * (q * 10 ** 3) * dio * lh) / (A_flow * g * lam)
                return xe
            elif hs == 2:
                xe = xi + (np.pi * (q2 * 10 **3 *doi) * lh) / (A_flow * g * lam)
                return xe
            else:
                if ch == 1:
                    xe = xi + (np.pi * ((q * 10 ** 3 * dio)+ (q2 * 10 ** 3 *doi)) * lh) / (A_flow * g * lam)
                else:
                    xe = xi + (np.pi * ((qav * 10 ** 3 * (dio + doi)) * lh)) / (A_flow * g * lam)
                return xe
        else:
            xe = xi + (qsol * C_heated * lh) / (C_flow * g * lam)
            return xe
    
    def cal_xi(self, q, q2, doi, dio, geo, hs, g, xe, lh, lam, ch = 1):
        # calculate heated section parameter
        if hs == 1:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * dio)
            R_heated = doi
            qsol = (q + q2) * (10 ** 3) # Lambda = kJ/kg
        elif hs == 2:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * doi)
            R_heated = 2 * doi
            qsol = (q + q2) * (10 ** 3)
        else:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * (dio + doi))        
            R_heated = 2 * doi
            qsol = (q + q2) * (10 ** 3)
        
        # calculate flow area
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        R_flow = doi * dio
        C_flow = (np.pi / 4) * doi ** 2
        
        # calculate quality
        if ch == 0:
            if geo == 'R':
                xi = xe - (qsol * R_heated * lh) / (R_flow * g * lam) 
                return xi
            elif geo == 'A':
                xi = xe - (np.pi*(q2 * (10 ** 3) * doi + q * (10 ** 3) * dio) * lh) / (A_flow * g * lam)
                return xi
            else:
                xi = xe - (qsol * C_heated * lh) / (C_flow * g * lam)
                return xi
        elif ch == 1:
            xi = xe- (qsol * C_heated * lh) / (C_flow * g * lam)
            return xi
        else:
            xi = xe- (q * (10 ** 3) * C_heated * lh) / (C_flow * g * lam)
            return xi

    def cal_ent(self, q, doi, dio, geo, hs, g, enthin, lh, lam):
        if hs == 1:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * dio)
            R_heated = doi
            qsol = q * (10 ** 3) # Lambda = kJ/kg
        elif hs == 2:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * doi)
            R_heated = 2 * doi
            qsol = (q) * (10 ** 3)
        else:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * (dio + doi))     
            R_heated = 2 * doi
            qsol = (q) * (10 ** 3)
        
        # calculate flow area
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        R_flow = doi * dio
        C_flow = (np.pi / 4) * doi ** 2

        if geo == 'R':
            enth = enthin + (qsol / R_flow) / (R_flow * g)
            return enth
        elif geo == 'A':
            enth = enthin + (qsol / A_flow) / (A_flow * g)
            return enth
        else:
            enth = enthin + (qsol / C_flow) / (C_flow * g)
            return enth


    def calXi1(self, cpf, dtin, lam): # Inlet thermal equilibrium quality
        Xi = - (cpf * dtin) / lam
        return Xi

    def calXi2(self, degsubin, lam):
        Xi = - (degsubin * 1e3) / lam
        return Xi

    def calGsat(self, q, doi, lh, cpf, dtin, dio, geo, hs, dh):
        qq = q * (10 ** 6)
        R_heated_1h = doi * lh
        R_heated_2h = 2 * doi * lh
        R_flow = doi * dio
        A_heated_1h = (np.pi * dio) * lh
        A_heated_2h = (np.pi * (dio + doi)) * lh
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        C_heated = (np.pi) * doi * lh
        C_flow = (np.pi / 4) * doi ** 2

        if geo == 'R':
            if hs == 1:
                Gsat = (qq * R_heated_1h) / (R_flow * cpf * dtin)
                return Gsat
            else:
                Gsat = (qq * R_heated_2h) / (R_flow * cpf * dtin)
                return Gsat
        elif geo == 'A':
            if hs == 1:
                Gsat = (qq * A_heated_1h) / (A_flow * cpf * dtin)
                return Gsat
            else:
                Gsat = (qq * A_heated_2h) / (A_flow * cpf * dtin)
                return Gsat
        else:
            Gsat = (qq * C_heated) / (C_flow * cpf * dtin)
            return Gsat

    def calDtOSV(self, q, g, cpf, dio, doi, lh, geo, hs, datap, tosv, tsat, ti):
        qq = q * (10 ** 6)
        R_heated_1h = doi * lh
        R_heated_2h = 2 * doi * lh
        R_flow = doi * dio
        A_heated_1h = (np.pi * dio) * lh
        A_heated_2h = (np.pi * (dio + doi)) * lh
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        C_heated = (np.pi) * doi * lh
        C_flow = (np.pi / 4) * doi ** 2

        if datap == 'OSV':
            dTOSV_Cal = tosv
        else:
            if geo == 'R':
                if hs == 1:
                    dTOSV_Cal = tsat - (ti + ((qq)*R_heated_1h) / (g*cpf*R_flow))
                    return dTOSV_Cal
                else:
                    dTOSV_Cal = tsat - (ti + ((qq)*R_heated_2h) / (g*cpf*R_flow))
                    return dTOSV_Cal
            elif geo == 'A':
                if hs == 1:
                    dTOSV_Cal = tsat - (ti + ((qq)*A_heated_1h) / (g*cpf*A_flow))
                    return dTOSV_Cal
                else:
                    dTOSV_Cal = tsat - (ti + ((qq)*A_heated_2h) / (g*cpf*A_flow))
                    return dTOSV_Cal
            else:
                dTOSV_Cal = tsat - (ti + ((qq)*C_heated / (g*cpf*C_flow)))
                return dTOSV_Cal

    def cal_de(self, geo, hs, doi, dio, dh, lh):
        if hs == 1:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * dio)
            R_heated = doi
        elif hs == 2:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * doi)
            R_heated = 2 * doi
        else:
            C_heated = (np.pi) * doi
            A_heated = (np.pi * (dio + doi))     
            R_heated = 2 * doi
        
        # calculate flow area
        A_flow = (np.pi / 4) * (doi ** 2 - dio ** 2)
        R_flow = doi * dio
        C_flow = (np.pi / 4) * doi ** 2

        if geo == 'R':
            De = 4 * R_flow / R_heated
            return round(De, 6), round(R_heated, 6), round(R_flow, 6)
        elif geo == 'A':
            De = 4 * A_flow / A_heated
            return round(De, 6), round(A_heated, 6), round(A_flow, 6)
        else:
            De = dh
            return round(De, 6), round(C_heated, 6), round(C_flow, 6)

    def calQ(self, doi, dio, lh, g, geo, hs, dh, ho, hi, de):
        q_guess = g * (de / (4 * lh)) * (ho - hi) * 1e-6

        # R_heated_1h = doi*lh
        # R_heated_2h = 2*doi*lh
        # R_flow = doi*dio
        # A_heated_1h = (np.pi*dio)*lh
        # A_heated_2h = (np.pi*(dio+doi))*lh
        # A_heated_3h = (np.pi*doi)
        # A_flow = (np.pi/4)*(doi**2-dio**2)
        # C_heated = (np.pi)*doi*lh
        # C_flow = (np.pi/4)*doi**2
        # R_phpw_1h = doi/(2*(doi+dio))
        # R_phpw_2h = doi/(doi+dio)
        # A_phpw_1h = dio/(2*(doi+dio))
        # A_phpw_3h = doi/(doi+dio)
        # param_DD = de/dh

        # if geo == 'R':
        #    if hs == 1:
        #        q_guess = g*param_DD*(R_flow/R_heated_1h)*(ho-hi)*1e-6
        #    else:
        #        q_guess = g*param_DD*(R_flow/R_heated_2h)*(ho-hi)*1e-6
        # elif geo == 'A':
        #    if hs == 1:
        #        q_guess = g*param_DD*(A_flow/A_heated_1h)*(ho-hi)*1e-6
        #    elif hs == 2:
        #        q_guess = g*(A_flow/A_heated_2h)*(ho-hi)*1e-6
        #    else:
        #        q_guess = g*param_DD*(A_flow/A_heated_3h)*(ho-hi)*1e-6
        # else:
        #    q_guess = g*(C_flow/C_heated)*(ho-hi)*1e-6
        return q_guess