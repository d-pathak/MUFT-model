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
mod_start = 60
beta = 1.55
Fr = 0.61 #derived in MUFT_flow_param.py
reaches = 3
#reach characteristics
L_list = [3130, 4660, 2990]
b_list = [0.1554, 0.0047, 0.0489]
c_list = [0.3967, 0.8699, 0.4352]

def calc_dQdt(y, t, Qi, reach):

    Lreach = L_list[reach]
    b = b_list[reach]
    c = c_list[reach]

    tf = int(mod_start + t)
    Tflow = ((Lreach / (b * (y**c))) / (m * (1 + beta))) / 60  #unit:min
    tau_fl = (1 - Fr) * Tflow #unit:min
    tlag = round(tf - int(tau_fl/tstep))
        
    dQdt = 5 * 60 * ((Qi[tlag] - y) / (Fr * Tflow * 60))  

    return dQdt
    

def unsteady_flow_routing(input_file):

    flow = pd.read_csv(input_file)
    flow['date_time'] =  pd.to_datetime(flow['date_time'], format='%d/%m/%Y %H:%M')

    Qi = flow['input_flow'].to_list()
    
    #store results and simulations
    routing = []
    final = []
    
    for reach in tqdm(range(reaches)):                 
                    
        if reach == 0:
            Qi = Qi
            y0 = Qi[mod_start]
        else:
            Qi = Qrouted
            y0 = Qi[mod_start]

        
        mod_len = len(Qi) - mod_start 
        t = np.arange(0,mod_len)

        routing = odeint(calc_dQdt, y0, t, args=(Qi,reach))
        Qrouted = np.array(routing[:,0]).tolist()
        Qout = np.array(routing[:,0]).tolist()

        final.append(Qout)
        Qout = []


    Qoutflow = pd.DataFrame(final).T
    final_df = Qoutflow
    final_df.to_csv('MUFT_Qout_routing.csv')
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