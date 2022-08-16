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
m = 5/3 #as per Chapra 2008
beta = 1.2 #solute-lag coefficient (estimate from conservative solute experiment)
tstep = 5 #model time-step, 5 min here
mod_start = int(24 * 60 / tstep) #to estimate from 6/8/2019 00:00 in the demo data, skip 1 day for t_lag
mod_len = int(24 * 60 / tstep) #for day-wise parameter estimation
krea = 0.0005 #reaeration coefficient
reaches = 2 #number of reaches in the network


def calc_dcdt(y, t, p, day, inputs, Ci): 

    try:        
        gppmax = p['gppmax'].value
        kpar = p['kpar'].value
        er = p['er'].value
    except:
        gppmax, kpar, er = p 

    tf = int(day * 288 + mod_start + t)
    t_lag = round(tf - int(inputs['tau_s'][tf]/tstep))

    gpp = (gppmax * inputs['PAR'][tf]) / (kpar + inputs['PAR'][tf]) 

    dCdt = tstep*(
        ((Ci[t_lag] - y) * (inputs['Qi'][t_lag] / (inputs['Qsite'][tf] * inputs['Tadz'][tf]))) + 
        (krea * (inputs['Cs'][tf] - y)) + 
        ((gpp - er) / inputs['dep'][tf])
        )
           
    return dCdt


#odeint, residual functions
def g(t, y0, p, day, inputs, Ci):
    soln = odeint(calc_dcdt, y0, t, args=(p, day, inputs, Ci))
    return soln[:,0]


def residual(p, t, y0, data, day, inputs, Ci):
    model = g(t, y0, p, day, inputs, Ci)
    return (model - data).ravel()


def solute_routing(input_file,networkfile):

    metab_inputs = pd.read_csv(input_file)
    network = pd.read_csv(networkfile)

    metab_inputs['date_time'] = pd.to_datetime(metab_inputs['date_time'],format='%d/%m/%Y %H:%M')
    datetimes = pd.to_datetime(metab_inputs['date_time'])
    metab_inputs['date'] = datetimes.dt.date
    metab_inputs.set_index(['reach'],inplace=True)

    metab_inputs['Catm'] = -0.00005858 * metab_inputs['Temp']**3 + 0.007195 * metab_inputs['Temp']**2 - 0.39509 * metab_inputs['Temp'] + 14.586
    metab_inputs['VP'] = 0.0000802 * metab_inputs['Temp']**3 - 0.000717 * metab_inputs['Temp']**2 + 0.0717 * metab_inputs['Temp'] + 0.539
    metab_inputs['Cs'] = (metab_inputs['Catm'] * (metab_inputs['P'] - metab_inputs['VP'])) / (101.325 - metab_inputs['VP'])

    reach_list = metab_inputs.index.unique()
    metab_inputs['L'] = np.nan
    metab_inputs['W'] = np.nan
    metab_inputs['a_dash'] = np.nan
    metab_inputs['b_dash'] = np.nan
    metab_inputs['tau_fl'] = np.nan

    for reach in reach_list:    
        metab_inputs['L'].mask(metab_inputs.index==reach,network['L'].iloc[reach-1],inplace=True) 
        metab_inputs['W'].mask(metab_inputs.index==reach,network['W'].iloc[reach-1],inplace=True)
        metab_inputs['a_dash'].mask(metab_inputs.index==reach,network['a_dash'].iloc[reach-1],inplace=True)
        metab_inputs['b_dash'].mask(metab_inputs.index==reach,network['b_dash'].iloc[reach-1],inplace=True)
        metab_inputs['tau_fl'].mask(metab_inputs.index==reach,network['tau_fl'].iloc[reach-1],inplace=True)

    metab_inputs['tau_s'] = m * (1 + beta) * metab_inputs['tau_fl'] / 60 #min
    metab_inputs['Tflow']= ((metab_inputs['a_dash'] * (metab_inputs['Qsite']**(metab_inputs['b_dash']-1))) + metab_inputs['tau_fl']) / 60 #obs/sim flow at reach end 
    metab_inputs['dep'] =  (metab_inputs['Qsite'] * (metab_inputs['Tflow']*60)) / (metab_inputs['W'] * metab_inputs['L']) 
    metab_inputs['Tsadz'] = (m * (1 + beta) * metab_inputs['Tflow'])
    metab_inputs['Tadz'] = metab_inputs['Tsadz'] - metab_inputs['tau_s']

    metab_inputs.reset_index(inplace=True)
    metab_inputs.set_index(['reach','date'],inplace=True)
    #date count
    tot_days = len(metab_inputs.index.get_level_values('date').unique()) - 1 # 1 date skip for t_lag
    date_list = pd.to_datetime(metab_inputs['date_time']).dt.date.unique().tolist()

    result = []
    store_results = []
    final = []
    sim = []

    for reach in tqdm(range(reaches)):

        filter = metab_inputs.loc[reach+1]
        filt1 = filter.reset_index(inplace=False)
        df = pd.DataFrame(filt1.iloc[-1]).T
        df1 = pd.DataFrame(np.repeat(df.values, 5, axis=0))
        df1.columns = df.columns
        df1.set_index(['date'],inplace=True)
        combine_df = [filter, df1]
        df_comb = pd.concat(combine_df)
        df_comb.reset_index(inplace=True)
        filt_inputs= df_comb[['Qi','Qsite','Ci','PAR','Cs','dep','Tadz','tau_s']]
        inputs = filt_inputs.to_dict()  
            
        if reach == 0:
            Ci = inputs['Ci']
        else:
            Ci = Cinarray


        for day in range(tot_days):

            date = date_list[day+reach+1]

            obs_filt = metab_inputs.loc[reach+1,date]
            data = obs_filt['Cobs']        

            if day == 0:
                y0 = Ci[mod_start]            
            else:
                y0 = solve2[-1]  
                             
            t = np.arange(0,mod_len)

            params = Parameters()
            params.add('gppmax', value = 0.006, min = 0) 
            params.add('kpar', min = 0.1, max = 0.5) 
            params.add('er', value = 0.002, min = 0) 

            #to specify reach-wise constraints
            # if reach == 0:
            #     params = Parameters()
            #     params.add('gppmax', value = 0.008, min = 0) 
            #     params.add('kpar', value = 0.1, min = 0.1, max = 0.5) 
            #     params.add('er', value = 0.005, min = 0)
            # else:
            #     params = Parameters()
            #     params.add('gppmax', value = 0.005, min = 0) 
            #     params.add('kpar', value = 0.2, min = 0.1, max = 0.5) 
            #     params.add('er', value = 0.001, min = 0)

            #nelder/leastsq method
            solve1 = minimize(residual, params, args=(t, y0, data, day, inputs, Ci), method='nelder')
            solve2 = data + solve1.residual.reshape(data.shape)

            #for rows in solve1:
            result.append(solve1)
            store_results.append(solve1)
            #for rows in solve2:
            final.append(solve2)
        
        
        final_array = np.asarray(final)
        Cinarray = final_array.flatten()

        arr = np.insert(Cinarray, 0, np.zeros((reach+1-1)*288) + np.nan) # adding NaN values at top to organise simulated data, compensate for 1 day lost at every reach.
        
        result = []
        final = []
        mc_result = []
        
        sim.append(arr)

        tot_days = tot_days - 1 #model run for n-1 days (n = number of simulation days for the preceding reach) to allow for t_lag

        #results
        stdoutOrigin=sys.stdout 
        sys.stdout = open('MUFTsolute_log_reach'+ str(reach+1) + '.txt', 'w')
        [report_fit(store_results[i]) for i in range(len(store_results))]
        sys.stdout.close()
        sys.stdout=stdoutOrigin

        store_results = []


    final_df = pd.DataFrame(sim).T
    final_df.to_csv('MUFT_generalsolute_sim.csv')

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
        solute_routing(f1,f2)


   
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print('Something went wrong {0}'.format(e))
