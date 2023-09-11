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
        #simdata['Day'] = simdata['Time'] % 365
        #simdata['Year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        #simdata['date'] = simdata.apply(
        #    lambda x: datetime.date(int(x['Year']), 1, 1) + datetime.timedelta(int(x['Day']) - 1), axis=1)

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
        sv = [s for s in self.sweep_variables]
        channels = [c for c in self.inset_channels]
        sv.append('Time')
        todo = {}
        for c in channels:
                todo[c] = 'mean'
        res = adf.groupby(sv, as_index=False).agg(todo)
        res.to_csv(os.path.join(self.working_dir, 'InsetChart.csv'), index=False)


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

class AgeDistributionAnalyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir='./',
                 burnin=None, filter_exists=False):

        super(AgeDistributionAnalyzer, self).__init__(working_dir=working_dir,
                                                    filenames=["output/MalariaSummaryReport_Daily.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.burnin = burnin
        self.filter_exists = filter_exists

    def filter(self, simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for fname in self.filenames:
            times = data[fname]['DataByTime']['Time Of Report']
            ages = data[fname]['Metadata']['Age Bins']
            d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin']
            time = []
            pop = []
            age = []
            for t in range(len(times)):
                for a in range(len(ages)):
                    time.append(times[t])
                    age.append(ages[a])
                    pop.append(d[t][a])
                    
            
            simdata = pd.DataFrame({'Time': time,
                                    'Age': age,
                                    'Pop': pop})
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

        adf.to_csv((os.path.join(self.working_dir, 'PopulationAgeStructure.csv')), index=False)


if __name__ == "__main__":

    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform

    
    expts = {#'exp_name' : 'exp_id'
             ######################
             #'sapone_simple_5_seed_30_degree' : '1433a53f-c830-48db-afc4-0eaae0afa75a',
             #'sapone_simple_5_seed_27_degree' : '1227520e-b5f7-40f6-ab18-433fc6720668',
             #'sapone_simple_5_seed_32_degree' : '8aa718ab-e366-48d7-bf9b-3feb65c2e6a7',
             #'sapone_simple_5_seed_32_degree_50yr' : 'b437ac26-d14f-42da-aa9f-bc9584085222',
             #'sapone_simple_5_seed_27_degree_50yr' : '654a1308-a4e6-4b3d-b7b4-cb47006524f9',
             #'sapone_simple_5_seed_30_degree_50yr' : '39fed015-4d5d-4c4c-85f6-7251eaa44273',
             #'test_calib_5yr': '4a4833ee-30d8-40ba-b2ce-870956da5ca6',
             'test_calib_50yr': '5c0f0421-443d-42cb-bb41-496a234703b6',
             #'test_calib_3yr': '6ea86674-0152-40bf-9354-ca9ddff7125b',
             
    }
    

    jdir =  '/projects/b1139/habitat_exploration/experiments'
    wdir_root =  '/projects/b1139/habitat_exploration/simulation_output'
    
    sweep_variables = ['Run_Number', 'x_Temporary_Larval_Habitat', 'scale_constant','scale_temp_rain','scale_water_veg'] 

    # set desired InsetChart channels to analyze and plot
    channels_inset_chart = ['Statistical Population', 'Adult Vectors', 'Daily Bites per Human', 'Daily EIR', 'Rainfall', 'Air Temperature', 'PCR Parasite Prevalence', 'New Clinical Cases']

    
    with Platform('SLURM_LOCAL',job_directory=jdir) as platform:

        for expname, exp_id in expts.items():
            wdir = os.path.join(wdir_root,expname)
            if not os.path.exists(wdir):
                os.mkdir(wdir)
            analyzer = [InsetChartAnalyzer(expt_name=expname,
                                      channels=channels_inset_chart,
                                      sweep_variables=sweep_variables,
                                      working_dir=wdir),
                        MonthlyAnalyzer(expt_name = expname,
                                        sweep_variables=sweep_variables,
                                        working_dir=wdir, 
                                        start_year=2020, 
                                        end_year=2023)]
            
            # Create AnalyzerManager with required parameters
            manager = AnalyzeManager(configuration={},ids=[(exp_id, ItemType.EXPERIMENT)],
                                     analyzers=analyzer, partial_analyze_ok=True)
            # Run analyze
            manager.analyze()
            


