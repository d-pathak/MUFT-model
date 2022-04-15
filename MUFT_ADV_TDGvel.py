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
from datetime import datetime
import calendar
import corner
import math
import argparse
import io
import errno


#parameter values, constants
m = 5/3 #as per Chapra 2008
mod_start = 163 #to estimate from 5/8/2019 00:00
tstep = 5 #model time-step
krea = 0.00025 #reaeration coefficient
days = 3 #period of estimation
start_day = 5 #first day number in the period of estimation 


def calc_dcdt(y, t, p, Qi, Qsite, PAR, Ci, Cs, dep, day, Tsadv, Tflowadv_emp):

    try:        
        gppmax = p['gppmax'].value
        kpar = p['kpar'].value
        er = p['er'].value
    except:
        gppmax, kpar, er = p 

    tf = int(day * 288 + mod_start + t) 
    alpha= Tsadv[tf] - Tflowadv_emp[tf]         
    t_lag = round(tf - int(alpha/tstep))
       
    gpp = (gppmax * PAR[tf]) / (kpar + PAR[tf]) 
    
    dCdt = tstep*(
        ((Ci[t_lag] - y) * (Qi[t_lag] / (Qsite[tf] * Tsadv[tf]))) + 
        (krea * (Cs[tf] - y)) + 
        ((gpp - er) / dep[tf])
        )            
    return dCdt


#odeint, residual functions
def g(t, y0, p, Qi, Qsite, PAR, Ci, Cs, dep, day, Tsadv, Tflowadv_emp):
    soln = odeint(calc_dcdt, y0, t, args=(p, Qi, Qsite, PAR, Ci, Cs, dep, day, Tsadv, Tflowadv_emp))
    return soln[:,0]


def residual(p, t, y0, data, Qi, Qsite, PAR, Ci, Cs, dep, day, Tsadv, Tflowadv_emp):
    model = g(t, y0, p, Qi, Qsite, PAR, Ci, Cs, dep, day, Tsadv, Tflowadv_emp)
    return (model - data).ravel()


def solute_model(metab_inputs):
    
    result = []
    store_results = []
    final = []
    sim = []
    
    
    Qi = metab_inputs['Qi'].tolist()
    Qsite = metab_inputs['Qsite'].tolist()
    Ci = metab_inputs['Ci'].tolist()
    PAR = metab_inputs['PAR'].tolist()
    Cs = metab_inputs['Cs'].tolist()
    dep = metab_inputs['dep'].tolist()
    temp = metab_inputs['Temp'].tolist()
    Tsadv = metab_inputs['Tsadv'].tolist()
    Tflowadv_emp = metab_inputs['Tflowadv_emp'].tolist()
    Tuadv = metab_inputs['Tuadv'].tolist()    

  
    for day in tqdm(range(days)):
    
        day_no = day + start_day
        
        obs_filt = metab_inputs.loc[1,day_no] #for 1 reach estimation
        data = obs_filt['Cobs'] 
                      
        mod_len = int(24 * 60 / tstep)                                
        t = np.arange(0,mod_len)
        if day == 0:
            y0 = Ci[mod_start]            
        else:
            y0 = solve2.values[-1] 

        #ADZ model: parameter constraints from the modified two-station model outputs
        params = Parameters()        
        params.add('gppmax', value = 0.0018, min = 0)
        params.add('kpar', min = 0.1, max = 0.5)
        params.add('er', value = 0.003, min = 0) 

        #nelder/leastsq method
        solve1 = minimize(residual, params, args=(t, y0, data, Qi, Qsite, PAR, Ci, Cs, dep, day, Tsadv, Tflowadv_emp), method='nelder')
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
    sys.stdout = open("MUFT_solute_log.txt", "w")
    report_fit(store_results[0])
    report_fit(store_results[1]) 
    report_fit(store_results[2])   
    sys.stdout.close()
    sys.stdout=stdoutOrigin

    dates_list = metab_inputs['date_time'].iloc[mod_start:mod_start+days*mod_len]
    PAR_list = metab_inputs['PAR'].iloc[mod_start:mod_start+days*mod_len]
    dates = pd.DataFrame(dates_list)
    PAR_values = pd.DataFrame(PAR_list)
    final_df['date_time'] =  dates['date_time'].values
    final_df['PAR'] =  PAR_values['PAR'].values
    final_df['date_time'] = pd.to_datetime(final_df['date_time'],format='%d/%m/%Y %H:%M')
    final_df=final_df.rename(columns={0:'DOsim','date_time':'dates','PAR':'PAR'})

    gpp_daily = np.full(len(final_df),np.nan)
    er_daily = np.full(len(final_df),np.nan)

    for day in range(days):

        gppmax = store_results[day].params['gppmax']
        kpar = store_results[day].params['kpar']
        er = store_results[day].params['er']
        
        for i in range(mod_len):
            gpp_daily[day*mod_len + i] = ((gppmax * final_df['PAR'][day*mod_len + i]) / (kpar + final_df['PAR'][day*mod_len + i])) * 60 * 24
            er_daily[day*mod_len + i] = er * 60 * 24


    final_df['gpp_daily'] = gpp_daily
    final_df['er_daily']= er_daily
    final_df['nep_daily'] = (final_df['gpp_daily'] - final_df['er_daily'])
    final_df.to_csv('MUFT_metab_out.csv',sep='\t')

    fig, ax = plt.subplots(1,figsize = (10,5))
    obs = metab_inputs['Cobs'].iloc[mod_start:mod_start+3*mod_len]
    ax.plot(final_df['dates'], final_df['DOsim'], c='blue')
    ax.scatter(final_df['dates'],obs,c='red',s=25)
    plt.savefig('MUFT_dofit.png')

    return final_df


def main(metab_inputs):
    solute_model(metab_inputs)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('inputsfile', help='Name of the input file you want to load')    
    args = parser.parse_args()

    with io.open(args.inputsfile, 'r', encoding='utf-8') as f1: 
        inputread = csv.reader(f1)
       
        metab_inputs = pd.read_csv(f1)
        L = 4660 #reach length
        W = 316 #reach width        
        metab_inputs['Catm'] = -0.00005858 * metab_inputs['Temp']**3 + 0.007195 * metab_inputs['Temp']**2 - 0.39509 * metab_inputs['Temp'] + 14.586
        metab_inputs['VP'] = 0.0000802 * metab_inputs['Temp']**3 - 0.000717 * metab_inputs['Temp']**2 + 0.0717 * metab_inputs['Temp'] + 0.539
        metab_inputs['Cs'] = (metab_inputs['Catm'] * (metab_inputs['P'] - metab_inputs['VP'])) / (101.325 - metab_inputs['VP'])
        #celerity, velocity calculations
        metab_inputs['us'] = 0.0047 * (metab_inputs['Qsite'] ** 0.8699)
        metab_inputs['uadv'] = metab_inputs['us'] 
        metab_inputs['cadv'] = m * metab_inputs['uadv']
        metab_inputs['Tflowadv_emp'] = (L / metab_inputs['cadv']) / 60 #unit:min
        metab_inputs['Tsadv'] = (L / metab_inputs['us']) / 60 #unit:min
        metab_inputs['Tuadv'] = (L / metab_inputs['uadv']) / 60 #unit:min
        metab_inputs['dep'] = (metab_inputs['Qsite'] * (metab_inputs['Tflowadv_emp']*60)) / (W * L)

        metab_inputs['date_time'] = pd.to_datetime(metab_inputs['date_time'],format='%d/%m/%Y %H:%M')
        datetimes = pd.to_datetime(metab_inputs['date_time'])
        metab_inputs['day'] = datetimes.dt.day
        metab_inputs.set_index(['reach','day'],inplace=True)
        metab_inputs.to_csv('MUFT_inputs.csv')
        main(metab_inputs)
        
