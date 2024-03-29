verbose = true
 if verbose && !isdefined(Main,:RvSpectMLBase)  println("# Loading RvSpecMLBase")    end
 using RvSpectMLBase
 if verbose && !isdefined(Main,:EchelleInstruments)  println("# Loading EchelleInstruments")    end
 using EchelleInstruments, EchelleInstruments.NEID
 using DataFrames, Query

# USER: You must create a data_paths.jl file in one of the default_paths_to_search listed below. It need only contain one line:
# solar_data_path = "/home/eford/Data/SolarSpectra/NEID_solar/"
target_subdir = "good_days/DRPv0.7/"   # USER: Replace with directory of your choice
 #fits_target_str = "Solar"
 fits_target_str = "Sun"
 output_dir = "examples/output"
 paths_to_search_for_param = [pwd(),"examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/group/ebf11/default/ebf11/neid_solar"]
 # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
 pipeline_plan = PipelinePlan()

reset_all_needs!(pipeline_plan)
if need_to(pipeline_plan,:read_spectra)
   if verbose println("# Finding what data files are avaliable.")  end
   eval(read_data_paths(paths_to_search=paths_to_search_for_param))
   @assert isdefined(Main,:solar_data_path)
   df_files = make_manifest(solar_data_path, target_subdir, NEID )
   #df_files = EchelleInstruments.make_manifest(solar_data_path, target_subdir, NEID )

   if verbose println("# Reading in customized parameters from param.jl.")  end
   eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))
   end

if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
   @time all_spectra = map(NEID.read_data,eachrow(df_files_use))
   dont_need_to!(pipeline_plan,:read_spectra)


#=
# Pre-ship corrections
if verbose println("# Applying wavelength corrections.")  end
   @assert isdefined(Main,:ancilary_solar_data_path)
   NEID.read_drift_corrections!(joinpath(ancilary_solar_data_path,"SolarRV20190918_JD_SciRV_CalRV.txt"), df_files_use)
   NEID.read_barycentric_corrections!(joinpath(ancilary_solar_data_path,"SolarTelescope2019-09-18_inclGravRedshiftAirMassAltitude.csv"), df_files_use)
   NEID.read_differential_extinctions!(joinpath(ancilary_solar_data_path,"20190918_diff_ex_full_fixed.txt"), df_files_use)
   apply_doppler_boost!(all_spectra,df_files_use)
   all_spectra
end
=#
