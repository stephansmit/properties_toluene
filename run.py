import fluids

if __name__=="__main__":
    Pref = 2e5
    Tref = 500
    gamma =1.06
    R = 90.23
    ig = fluids.IdealGasFluid("Toluene", Pref, Tref)
    ig.create_table(1e5, 2e5, 400, 500, 10, 10)

    print(ig.df)

