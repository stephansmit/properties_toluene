import CoolProp.CoolProp as CP
from fluidmodels_su2.SU2Models import *
import numpy as np
import pandas as pd
class Fluid(object):
    def __init__(self, name):
        self.name = name
    
    def create_table(self,Pmax, Pmin,Tmax, Tmin, Np, Nt):
         P_vec= np.geomspace(Pmin, Pmax, Np)
         T_vec = np.linspace(Tmin, Tmax, Nt)
         T = np.array([[t for t in T_vec ] for p in P_vec])
         P = np.array([[p for t in T_vec ] for p in P_vec])
         self.df = pd.DataFrame(data= [T.flatten(), P.flatten()]).transpose()
         self.df.columns = ["T", "P"]
         self.df['R']= self.df.apply(lambda x: self.get_density(x['P'], x['T']),axis=1)
         self.df['S']= self.df.apply(lambda x: self.get_entropy(x['P'], x['T']),axis=1)
         self.df['E']= self.df.apply(lambda x: self.get_internal_energy(x['P'], x['T']),axis=1)


class SU2Fluid(Fluid):
    def __init__(self, name):
        Fluid.__init__(self,name)

    def get_thermal_conductivity(self, P, T):
        self.Model.SetTDState_PT(P,T)
        return self.Model.GetThermalConductivity()
    
    def get_soundspeed(self, PT):
        self.Model.SetTDState_PT(P,T)
        return self.Model.GetSoundSpeed()

    def get_density(self, P, T):
        self.Model.SetTDState_PT(P,T)
        return self.Model.GetDensity()

    def get_internal_energy(self, P, T):
        self.Model.SetTDState_PT(P,T)
        return self.Model.GetStaticEnergy()

    def get_entropy(self, P, T):
        self.Model.SetTDState_PT(P,T)
        return self.Model.GetEntropy()

    def get_viscosity(self, P, T):
        self.Model.SetTDState_PT(P,T)
        return self.Model.GetLaminarViscosity()

class CoolPropFluid(Fluid):
    def __init__(self, name):
        Fluid.__init__(self,name)

    def get_thermal_conductivity(self, P, T):
        return CP.PropsSI("L", "P", P, "T", self.name)

    def get_density(self, P, T):
        return CP.PropsSI("D", "P", P, "T", self.name)

    def get_enthalpy(self, P, T):
        return CP.PropsSI("H", "P", P, "T", self.name)

    def get_entropy(self, P, T):
        return CP.PropsSI("D", "P", P, "T", self.name)

    def get_viscosity(self, P, T):
        return CP.PropsSI("V", "P", P, "T", self.name)


class IdealGasFluid(SU2Fluid):
    def __init__(self, name, Pref, Tref):
        SU2Fluid.__init__(self,name)
        self.R = CP.PropsSI("GAS_CONSTANT", self.name)* \
                 (CP.PropsSI("DMOLAR", "P", Pref, "T", Tref, "Toluene")/CP.PropsSI("D", "P", Pref, "T", Tref, self.name))
        self.gamma = CP.PropsSI("CPMASS", "P", Pref, "T", Tref, self.name)/CP.PropsSI("CVMASS", "P", Pref, "T", Tref, self.name)
        self.set_model()

    def set_model(self):
        self.Model = CIdealGas(self.gamma, self.R)

class SpecificIdealGasFluid(SU2Fluid):
    def __init__(self, name, gamma, R):
        SU2Fluid.__init__(self,name)
        self.R = R
        self.gamma = gamma
        self.set_model()

    def set_model(self):
        self.Model = CIdealGas(self.gamma, self.R)

class PengRobinsonFluid(IdealGasFluid):
    def __init__(self, name, Pref, Tref):
        self.Tcrit = CP.PropsSI("Tcrit", name)
        self.Pcrit = CP.PropsSI("Pcrit", name)
        self.acentric = CP.PropsSI("ACENTRIC", name)
        IdealGasFluid.__init__(self,name, Pref, Tref)

    def set_model(self):
        self.Model = CPengRobinson(self.gamma, self.R, self.Tcrit, self.Pcrit, self.acentric)

class SpecificPengRobinsonFluid(SpecificIdealGasFluid):
    def __init__(self, name, gamma, R):
        self.Tcrit = CP.PropsSI("Tcrit", name)
        self.Pcrit = CP.PropsSI("Pcrit", name)
        self.acentric = CP.PropsSI("ACENTRIC", name)
        SpecificIdealGasFluid.__init__(self,name, gamma, R)

    def set_model(self):
        self.Model = CPengRobinson(self.gamma, self.R, self.Tcrit, self.Pcrit, self.acentric)


