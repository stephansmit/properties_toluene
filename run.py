import fluids
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import re
plt.rc('text', usetex=True)

def get_columns(name, columns):
    return filter(lambda x: re.match('^'+name, x) is not None, columns)

def plot_figure(grouped, columnname, ylabel, title, xlim, ylim, figure_name, styles, logy, show_title):
    fig,ax = plt.subplots()
    for name,  group in grouped:
        style = styles[float(name)]
        grouped.get_group(name).plot(x='T', y=columnname+'_CP',ax = ax, color='r', linestyle=style, logy=logy)
        grouped.get_group(name).plot(x='T', y=columnname+'_PR',ax = ax, color='b', linestyle=style, logy=logy)
        grouped.get_group(name).plot(x='T', y=columnname+'_IG',ax = ax, color='g', linestyle=style, logy=logy)
    
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlabel(r'$\displaystyle T$ [K]', fontsize=16)
    ax.tick_params(labelsize=14)
    if show_title:
        ax.set_title(title, fontsize=20)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    legend1 = plt.legend(['CoolProp', 'PR', 'IG'], loc=4)
    ax.add_artist(legend1)
    handles = [mlines.Line2D([], 
                             [], 
                             color='k',  
                             linestyle=value,
                             markersize=10, 
                             label='P='+str(name)) 
               for name, value in styles.items()]
    handles_sorted = []
    handles_sorted.append(handles[2])
    handles_sorted.append(handles[1])
    handles_sorted.append(handles[0])
    handles_sorted.append(handles[3])
    legend2= plt.legend(handles=handles_sorted, loc=3)
    fig.savefig(figure_name)
    plt.close()

if __name__=="__main__":
    Pref = 2e5
    Tref = 500
    gamma =1.06
    R = 90.23



    ig = fluids.IdealGasFluid("Toluene", Pref, Tref)
    pr = fluids.PengRobinsonFluid("Toluene", Pref, Tref)
    cp = fluids.CoolPropFluid("Toluene")
    
    
    Pvec = np.array([0.1, 5, 10, 15])
    styles = {0.1:'-',5: '-.',10: '--', 15:':'}
    ig.create_table(Pvec, 500, 600, 10 )
    pr.create_table(Pvec, 500, 600, 10 )
    cp.create_table(Pvec, 500, 600, 10 )

    pr.df.columns = [str(col) + '_PR' for col in pr.df.columns]
    ig.df.columns = [str(col) + '_IG' for col in ig.df.columns]
    cp.df.columns = [str(col) + '_CP' for col in cp.df.columns]
    df = pd.concat([pr.df, ig.df, cp.df], axis=1)
    df['P']=df['P_PR'].apply(lambda x: "{0:0.1f}".format(x*1e-5))
    df['T']=df['T_PR']
#    for col in get_columns('k',df.columns):        df[col] = df[col]*1e3;
#    for col in get_columns('mu',df.columns):       df[col] = df[col]*1e6;
    for col in get_columns('dSdP_R',df.columns): df[col] = df[col]*1e3;
    for col in get_columns('dSdR_P',df.columns): df[col] = df[col]/-1e3;
    for col in get_columns('dHdR_P',df.columns): df[col] = df[col]/-1e3;



    grouped = df.groupby("P")
    plot_figure(grouped, 'R', r'$\displaystyle \rho$ [kg/m^3]',"Density with different models", [500,600], [0,1000], "density.png", styles,True, False)
    plot_figure(grouped, 'A', r'$\displaystyle c$ [m/s]',"Speed of sound with different models", [500,600], [150,240], "speedofsound.png", styles,False, False)
#    plot_figure(grouped, 'L', 'k [kW /m K]', "Thermal conductivity with different models", [500,600], [25,45], "thermalconductivity.png", plist, styles, False)
#    plot_figure(grouped, 'V', r'$\displaystyle \mu$ [$\displaystyle \mu$ Pa/ s]',"Viscosity with different models", [500,600], [10,15], "viscosity.png", plist, styles,False)
    plot_figure(grouped, 'dSdR_P', r'$\displaystyle -\frac{\partial s}{\partial \rho}|_P$ ',#[J m^3/kg^2 K]',
            "Parial derivative of entropy with respect to density at constant pressure", [500,600], [0.01,13], "dsdrho_p.png", styles,True, False)
    plot_figure(grouped, 'dSdP_R', r'$\displaystyle \frac{\partial s}{\partial P}|_\rho$',# [MJ/kg K Pa]',
            "Parial derivative of entropy with respect to pressure at constant density", [500,600], [0.5,300],"dsdp_rho.png", styles,True, False)
    plot_figure(grouped, 'dHdR_P', r'$\displaystyle -\frac{\partial h}{\partial \rho}|_P$',# [J m^6/kg^2]',
            "Parial derivative of enthalpy with respect to density at constant pressure", [500,600], [1,7500], "dhdrho_p.png", styles,True, False)
    plot_figure(grouped, 'dHdP_R', r'$\displaystyle \frac{\partial h}{\partial P}|_\rho$',# [KJ/kg Pa]',
            "Parial derivative of enthalpy with respect to pressure at constant density", [500,600], [0.2,130],"dhdp_rho.png", styles,True, False)

