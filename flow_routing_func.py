
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters, Parameter, report_fit, Minimizer, conf_interval, printfuncs
from scipy.integrate import odeint
from scipy import interpolate
import csv
import numdifftools
from pylab import shape
from tqdm import tqdm
import emcee
import datetime
import corner
import argparse
import io


def calc_dQdt(y, t, p, Qi):

    try:
        a_dash = p['a_dash'].value
        b_dash = p['b_dash'].value
        tfl = p['tfl'].value        
    except:
        a_dash, b_dash, tfl = p
    
    tf = int(t)

    #3 param celerity model, n= 1 is assumed
    tflow = ((a_dash * (y**(b_dash-1))) + tfl)
    dQdt = 3600 * ((Qi[tf] - y) / tflow)  

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
    L_reach = flow['reach_length'].to_list()

    #store results and simulations
    result = []
    store_results = []
    final = []
    sim = []
    
    mod_len = len(Qi) - 1
    data = flow['observed_flow'][:mod_len]
    t = np.arange(0,mod_len)
    y0 = Qi[0]

    params = Parameters()
    params.add('a_dash', value= 10000, min = 0) 
    params.add('b_dash', value = 0.5, min = 0, max = 1) 
    params.add('tfl', value = 10800, min = 0) 

    #original nelder/leastsq method
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
    ts = np.linspace(0,192,193)
    ax.plot(ts, final_df[0], c='blue')
    ax.scatter(ts,data,c='red',s=25)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Flow ($m^3$ $s^{-1}$)')
    plt.savefig('MUFT_flow__fit.png')
    from sklearn.metrics import r2_score 
    r2 = r2_score(final_df[0],data)
    r2

    sim = final_df.rename(columns = {0:'simulated_flow'})
    sim['date_time'] = flow['date_time']
    sim['date_time'] =  pd.to_datetime(sim['date_time'], format='%d/%m/%Y %H:%M')
    series = sim.set_index(sim['date_time'])
    del series['date_time']
    upsampled = series.resample('5Min')
    interpolated = upsampled.interpolate(method='linear')
    interpolated.to_csv('MUFT_simflow.csv', sep='\t')

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