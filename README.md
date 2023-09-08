# FE Project Template

Basic instructions:

1. Create a folder within simulation_inputs for a new site (ex. "sapone")
2. Create a new demographics.csv for your site(s)/node(s)
3. Run site_setup with step = 'setup_climate' to create a folder with climate files specifying:
     - date range to request from COMPS
     - Constant temperature value (optional)
     - rainfall shift (optional)
5. Populate other input files with information relevant to your site simulations:
     - vectors.csv - for vector species mix and parameters

       If applicable:
     - interventions_CM.csv - For treatment-seeking / case management
     - interventions_SMC.csv - For seasonal chemoproprevention 
     - interventions_IRS.csv - For indoor residual spraying
     - interventions_ITN.csv - For usage-dependent bednet distributions, also include:
         - interventions_ITN_age.csv     - age pattern of net usage
         - interventions_ITN_season.csv  - seasonal pattern of net usage
6. Run site_setup with step = 'burnin'. Copy experiment-id to burnin_id
7. Run site_setup with step = 'pickup'. Copy experiment-id to pickup_id
8. Run site_setup with step = 'calibrate' or 'analyze_pickup'. Results will be stored in simulation_outputs/site_name/exp_label
9. **For now** use plot_site_setup.Rmd to visualize outputs
