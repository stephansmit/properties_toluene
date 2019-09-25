import fluids

if __name__=="__main__":
    Pref = 2e5
    Tref = 500
    gamma =1.06
    R = 90.23
    ig = fluids.IdealGasFluid("Toluene", Pref, Tref)
    sig = fluids.SpecificIdealGasFluid("Toluene", gamma, R)
    spr = fluids.SpecificPengRobinsonFluid("Toluene", gamma, R)
    print(ig.gamma)
    print(ig.R)
    print(sig.gamma)
    print(sig.R)
    print(spr.acentric)

