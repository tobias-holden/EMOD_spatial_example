import os
import datetime
import pandas as pd
import numpy as np
import sys
import re
import random
from idmtools.entities import IAnalyzer	
from idmtools.entities.simulation import Simulation

## For plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
    

class InsetChartAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, sweep_variables=None, channels=None, working_dir=".", start_year=1970, years_to_keep=int(2)):
        super(InsetChartAnalyzer, self).__init__(working_dir=working_dir, filenames=["output/InsetChart.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.inset_channels = channels or ['Statistical Population', 'Adult Vectors', 'Daily Bites per Human', 'Daily EIR', 'Rainfall', 'Air Temperature']
        self.expt_name = expt_name
        self.start_year = start_year
        self.years_to_keep = years_to_keep

    def map(self, data, simulation: Simulation):
        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.inset_channels})
        simdata['Time'] = simdata.index
        cutoff = np.max(simdata['Time'])-self.years_to_keep*365
        simdata = simdata[simdata['Time']>=cutoff]
        simdata['Day'] = simdata['Time'] % 365
        simdata['Year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        simdata['date'] = simdata.apply(
            lambda x: datetime.date(int(x['Year']), 1, 1) + datetime.timedelta(int(x['Day']) - 1), axis=1)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            elif sweep_var == 'Run_Number' :
                simdata[sweep_var] = 0
        return simdata

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir)):
            os.mkdir(os.path.join(self.working_dir))

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(os.path.join(self.working_dir, 'InsetChart.csv'), index=False)


class MonthlyAnalyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir='./',
                 burnin=None, filter_exists=False, years_to_keep = int(2)):

        super(MonthlyAnalyzer, self).__init__(working_dir=working_dir,
                                                    filenames=["output/MalariaSummaryReport_Monthly_AllAge.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.burnin = burnin
        self.filter_exists = filter_exists
        self.years_to_keep = years_to_keep

    def filter(self, simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for fname in self.filenames:
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:(self.years_to_keep*12)]
            pfpr = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:(self.years_to_keep*12)]
            clinical_cases = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:(self.years_to_keep*12)]
            severe_cases = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:(self.years_to_keep*12)]
            pop = [x[1] for x in d]
            d = data[fname]['DataByTime']['PfPR_2to10'][:(self.years_to_keep*12)]
            PfPR_2to10 = d
            d = data[fname]['DataByTime']['Annual EIR'][:(self.years_to_keep*12)]
            annualeir = d
            simdata = pd.DataFrame({'month': range(1, (self.years_to_keep*12 + 1)),
                                    'PfPR': pfpr,
                                    'Cases': clinical_cases,
                                    'Severe_cases': severe_cases,
                                    'Pop': pop,
                                    'PfPR_2to10': PfPR_2to10,
                                    'annualeir': annualeir})
            adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                try:
                    adf[sweep_var] = simulation.tags[sweep_var]
                except:
                    adf[sweep_var] = '-'.join([str(x) for x in simulation.tags[sweep_var]])
            elif sweep_var == 'Run_Number':
                adf[sweep_var] = 0

        return adf

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return



        print(f'\nSaving outputs to: {os.path.join(self.working_dir)}')

        adf = pd.concat(selected).reset_index(drop=True)

        adf.to_csv((os.path.join(self.working_dir, 'AllAge_MonthlySummaryReport.csv')), index=False)
