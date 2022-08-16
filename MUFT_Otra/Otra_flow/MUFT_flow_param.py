
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters, Parameter, report_fit, Minimizer, conf_interval, printfuncs
from scipy.integrate import odeint
import csv
import numdifftools
from pylab import shape
from tqdm import tqdm
import argparse
import io

#parameter values, constants
m = 5/3
tstep = 5
mod_start = 288 #to estimate from 4/8/2019 00:00
beta = 1.55 #estimated using TDG vel
#reach characteristics
L = [3130, 4660, 2990]
b = [0.1554, 0.0047, 0.0489]
c = [0.3967, 0.8699, 0.4352]

def calc_dQdt(y, t, p, Qi):

    try:
        Fr = p['Fr'].value         
    except:
        Fr = p 
    
    tf = int(mod_start + t)

    Ts1 = L[0] / (b[0] * (y**c[0])) #unit: seconds
    Ts2 = L[1] / (b[1] * (y**c[1])) 
    Ts3 = L[2] / (b[2] * (y**c[2]))
    Ts = (Ts1 + Ts2 + Ts3) / 60   #unit:min 

    Tflow = Ts / (m * (1 + beta)) #unit:min    
    tau_fl = (1 - Fr) * Tflow  #unit:min
    tlag = round(tf - int(tau_fl/tstep))

    dQdt = 5 * 60 * ((Qi[tlag] - y) / (Fr * Tflow * 60)) 

    return dQdt


#odeint, residual functions
def g(t, y0, p, Qi):
    soln = odeint(calc_dQdt, y0, t, args=(p,Qi))
    return soln[:,0]


def residual(p, ts, y0, data, Qi):
    model = g(ts, y0, p, Qi)
    return (model - data).ravel()


def unsteady_flow_routing(input_file):

    flow = pd.read_csv(input_file)
    flow['date_time'] =  pd.to_datetime(flow['date_time'], format='%d/%m/%Y %H:%M')

    Qi = flow['input_flow'].to_list()
    
    #store results and simulations
    result = []
    store_results = []
    final = []
    sim = []
    
    mod_len = len(Qi) - 1 - mod_start
    data = flow['observed_flow'][mod_start:mod_start + mod_len]
    t = np.arange(0,mod_len)
    y0 = Qi[mod_start]

    params = Parameters()
    params.add('Fr', value = 0.6, min = 0, max = 1) 
    
    solve1 = minimize(residual, params, args=(t, y0, data, Qi), method='nelder')
    solve2 = data + solve1.residual.reshape(data.shape)

    #for rows in solve1:
    result.append(solve1)
    store_results.append(solve1)
    #for rows in solve2:
    final.append(solve2)
    final_array = np.asarray(final)
    Cinarray = final_array.flatten()
    sim.append(Cinarray)

    final_df = pd.DataFrame(sim).T

    stdoutOrigin=sys.stdout 
    sys.stdout = open("MUFT_flow_log.txt", "w")
    report_fit(store_results[0])
    sys.stdout.close()
    sys.stdout=stdoutOrigin

    fig, ax = plt.subplots(1,figsize = (10,5))
    ts = np.linspace(0,mod_len - 1, mod_len)
    ax.plot(ts, final_df[0], c='blue')
    ax.scatter(ts,data,c='red',s=25)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Flow ($m^3$ $s^{-1}$)')      
    plt.savefig('MUFT_flow__fit.png')

    final_df.to_csv('MUFT_simflow_Hekni.csv')
    
    return final_df


def main():
    # Create the parser
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Name of the file you want to load')
    args = parser.parse_args()

    with io.open(args.filename, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        unsteady_flow_routing(f)

   
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print('Something went wrong {0}'.format(e))