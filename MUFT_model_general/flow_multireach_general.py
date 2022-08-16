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


def calc_dQdt(y, t, p, Qi):

    try:
        a_dash = p['a_dash'].value
        b_dash = p['b_dash'].value
        tau_fl = p['tau_fl'].value        
    except:
        a_dash, b_dash, tau_fl = p
    
    tf = int(t)

    #3 param celerity model, n= 1 is assumed
    tflow = ((a_dash * (y**(b_dash-1))) + tau_fl)
    dQdt = 60 * 60 * ((Qi[tf] - y) / tflow)  #for hourly scale flow data

    return dQdt


#odeint, residual functions
def g(t, y0, p, Qi):
    soln = odeint(calc_dQdt, y0, t, args=(p,Qi))
    return soln[:,0]


def residual(p, ts, y0, data, Qi):
    model = g(ts, y0, p, Qi)
    return (model - data).ravel()


def unsteady_flow_routing(input_file,networkfile):

    inputdata = pd.read_csv(input_file)
    inputdata['date_time'] =  pd.to_datetime(inputdata['date_time'], format='%d/%m/%Y %H:%M')
    inputdata.set_index(['reach'],inplace=True)

    Qi = inputdata['input_flow'].to_list()
    
    network = pd.read_csv(networkfile)
    L = network['L'].tolist()

    #store results and simulations
    result = []
    store_results = []
    final = []
    sim = []
    simflow = []

    for reach in tqdm(range(len(L))):

        reach = reach + 1

        reach_data = inputdata.loc[reach]
        Qi_input = reach_data['input_flow'].to_list()
        data = reach_data['observed_flow']

        if reach == 1:
            Qi = Qi_input
            y0 = Qi[0]
        else:
            Qi = simflow
            y0 = simflow[0]
    
        mod_len = len(Qi)
        t = np.arange(0,mod_len)
        
        if reach == 1:
            params = Parameters()
            params.add('a_dash', value= 8000, min = 0) 
            params.add('b_dash', value = 0.4, min = 0.1, max = 1) 
            params.add('tau_fl', value = 3000, min = 0)
        else:
            params = Parameters()
            params.add('a_dash', value= 40000, min = 0) 
            params.add('b_dash', value = 0.5, min = 0.1, max = 1) 
            params.add('tau_fl', value = 8000, min = 0) 

        #nelder/leastsq method
        solve1 = minimize(residual, params, args=(t, y0, data, Qi), method='nelder')
        solve2 = data + solve1.residual.reshape(data.shape)

        simflow = solve2.tolist()

        #for rows in solve1:
        result.append(solve1)
        store_results.append(solve1)
        #for rows in solve2:
        final.append(solve2)
        final_array = np.asarray(final)
        Cinarray = final_array.flatten()
        final = []
        sim.append(Cinarray)

        stdoutOrigin=sys.stdout 
        sys.stdout = open('MUFTflow_log_reach' + str(reach) + '.txt', 'w')
        report_fit(store_results[reach-1])
        sys.stdout.close()
        sys.stdout=stdoutOrigin


    final_df = pd.DataFrame(sim).T
    final_df = final_df.rename(columns = {0:'reach1_Qsim',1:'reach2_Qsim'})
    reach1_df = inputdata.loc[1]
    reach1_df = reach1_df.reset_index()
    final_df['date_time'] = reach1_df['date_time']
    final_df['input_flow'] = reach1_df['input_flow']
    final_df['date_time'] =  pd.to_datetime(final_df['date_time'], format='%d/%m/%Y %H:%M')
    series = final_df.set_index(final_df['date_time'])
    del series['date_time']
    upsampled = series.resample('5Min')
    interpolated = upsampled.interpolate(method='linear')
    interpolated.to_csv('MUFT_generalflow_sim.csv', sep='\t')
    
    return final_df


def main():
    # Create the parser
    parser = argparse.ArgumentParser()
    parser.add_argument('inputdata', help='Name of the input data file')
    parser.add_argument('network', help='Name of the network description file')
    args = parser.parse_args()         
        

    with io.open(args.inputdata, 'r', encoding='utf-8') as f1, io.open(args.network, 'r', encoding='utf-8') as f2:
        reader = csv.reader(f1)
        network = csv.reader(f2)
        unsteady_flow_routing(f1,f2)


   
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print('Something went wrong {0}'.format(e))