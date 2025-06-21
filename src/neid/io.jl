"""
   IO functions for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

"""Create Dataframe containing filenames and key data for all files neid*.fits in directory"""
function make_manifest(data_path::String ; max_spectra_to_use::Int = 1000 )
    df_filenames = EchelleInstruments.make_manifest(data_path,r"^neidL[1,2]_\d+[T\.]\d+\.fits$")
    #=
    dir_filelist = readdir(data_path,join=true)
    idx_spectra = map(fn->occursin(r"^neidL1_\d+[T\.]\d+\.fits$", last(split(fn,'/')) ),dir_filelist)
    spectra_filelist = dir_filelist[idx_spectra]
    =#
    #=
    df_files = DataFrame(Filename = String[], target = String[], bjd = Float64[], ssbz=Float64[] )
    map(fn->add_metadata_from_fits!(df_files,fn),spectra_filelist)
    df_files
    =#
    #@assert length(spectra_filelist) >= 1
    #df_files = DataFrame(read_metadata(spectra_filelist[1]))
    df_files = DataFrame(read_metadata(df_filenames.Filename[1]))
    keys = propertynames(df_files)
    allowmissing!(df_files, keys[map(k->k∉[:Filename, :bjd, :target, :airmass],keys)] )
    #println("df_files = ", df_files)


    #if length(spectra_filelist) >= 2
    if length(df_filenames.Filename) >= 2
        max_idx = min(max_spectra_to_use,length(df_filenames.Filename))
        map(fn->add_metadata_from_fits!(df_files,fn),df_filenames.Filename[2:max_idx])
    end
    #=
    for i in 2:length(spectra_filelist)
        try
            add_metadata_from_fits!(df_files,spectra_filelist[i])

        catch
            println("# Problem reading ", spectra_filelist[i])
            fields_to_save = metadata_symbols_default(NEID2D())
            fields_str_to_save = metadata_strings_default(NEID2D())

            metadata_tmp = read_metadata_from_fits(spectra_filelist[i],fields=fields_to_save,fields_str=fields_str_to_save)
            println("Before: ", metadata_tmp)
            #println("After : ", metadata_tmp)
            #push!(df_files,metadata_tmp)
        end
    end
    =#
    df_files

end

"""Create Dict containing filename and default metadata from file."""
function read_metadata(fn::String)
    f = FITS(fn)
    read_metadata(f)
end

function read_metadata(f::FITS)
    fields_to_save = metadata_symbols_default(NEID2D())
    fields_str_to_save = metadata_strings_default(NEID2D())
    dict = read_metadata_from_fits(f,fields=fields_to_save,fields_str=fields_str_to_save)
    check_metadata_fields_expected_present(dict,fields_to_save)
    dict[:Filename] = f.filename
    if haskey(dict,:airmass) && typeof(dict[:airmass]) == String
       dict[:airmass] = parse(Float64,dict[:airmass])
    end
    hdr1 = read_header(f[1])
    if hdr1["DATALVL"] == 2
        dict_ccf = read_metadata_ccf(f)
        dict_expmeter = get_exposure_meter_summary(f)
        dict_activity = read_activity_table(f)
        dict = merge(dict,dict_ccf,dict_expmeter,dict_activity)
    end
    return dict
end

function read_metadata_ccf(f::FITS)
    fields_to_save = [ :drp_extsnr, :drp_ccfjdmod, :drp_ccfrvmod, :drp_dvrmsmod]
    #fields_str_to_save = [ "EXTSNR", "CCFJDMOD", "CCFRVMOD", "DVRMSMOD"]
    fields_str_to_save = [ "EXTSNR", "CCFJDMOD", "CCFRVMOD", "DVRMSMOD","BISMOD","EBISMOD","FWHMMOD"]
    dict = read_metadata_from_fits(f,fields=fields_to_save,fields_str=fields_str_to_save, hdu="CCFS")
end

function read_metadata_ccf(fn::String)
    f = FITS(fn)
    read_metadata_ccf(f)
end

function read_activity_table(fn::String)
    f = FITS(fn)
    read_activity_table(f)
end

"""Create Dict containing filename and default metadata from file."""
function read_activity_table(f::FITS)
    activity_index = read(f["ACTIVITY"],"INDEX")
    activity_value = read(f["ACTIVITY"],"VALUE")
    activity_uncert = read(f["ACTIVITY"],"UNCERTAINTY")
    telluric_zenith = read_header(f["TELLURIC"])["ZENITH"]
    telluric_wvapor = read_header(f["TELLURIC"])["WVAPOR"]
    Dict(vcat(Symbol.(activity_index) .=> read(f["ACTIVITY"],"VALUE"),
              Symbol.("σ_" .* activity_index) .=> read(f["ACTIVITY"],"UNCERTAINTY"),
              Symbol("telluric_zenith") => telluric_zenith,
              Symbol("telluric_wvapor") => telluric_wvapor ) )
end

"""Create Dict containing filename and default metadata from file."""
function read_exposure_meter(f::FITS)
    hdr1 = read_header(f[1])
    exptime = hdr1["EXPTIME"]
    λ_range = 38:97
    fiber = 5
    #expmeter_data = read(f["EXPMETER"],t_range,λ_range,fiber)
    expmeter_data = read(f["EXPMETER"]) # ,:,λ_range,fiber)
    t_len = size(expmeter_data,1)
    t_range = 2:(min(round(Int64,exptime)-1,t_len))
    if t_len < floor(Int64,exptime)
       @warn("*** Exposure meter contains fewer readings than floor(EXPTIME) in " * string(f.filename) )
    end
    view(expmeter_data, t_range,λ_range,fiber)
end

function read_exposure_meter(fn::String)
    f = FITS(fn)
    read_exposure_meter(f)
end

function get_exposure_meter_summary(f::Union{FITS,String})
    try
       expmeter_data = read_exposure_meter(f)
       mean_expmeter = mean(sum(expmeter_data,dims=2),dims=1)[1,1]
       mean_expmeter_blue = mean(sum(view(expmeter_data,:,1:20),dims=2),dims=1)[1,1]
       mean_expmeter_green = mean(sum(view(expmeter_data,:,21:40),dims=2),dims=1)[1,1]
       mean_expmeter_red = mean(sum(view(expmeter_data,:,41:60),dims=2),dims=1)[1,1]
       rms_expmeter = sqrt(var(sum(expmeter_data,dims=2),mean=mean_expmeter,corrected=false))
       rms_expmeter_blue = sqrt(var(sum(view(expmeter_data,:,1:20),dims=2),mean=mean_expmeter,corrected=false))
       rms_expmeter_green = sqrt(var(sum(view(expmeter_data,:,21:40),dims=2),mean=mean_expmeter,corrected=false))
       rms_expmeter_red = sqrt(var(sum(view(expmeter_data,:,41:60),dims=2),mean=mean_expmeter,corrected=false))
       
       expmeter_data_winsor = similar(expmeter_data)
       for i in 1:size(expmeter_data_winsor,2)
           expmeter_data_winsor[:,i] .= winsor(view(expmeter_data,:,i),count=2)
       end
       expmeter_mean_winsor = mean(sum(expmeter_data_winsor,dims=2),dims=1)[1,1]
       expmeter_rms_winsor = sqrt(var(sum(expmeter_data_winsor,dims=2),mean=expmeter_mean_winsor,corrected=false))
       return Dict(:expmeter_mean => mean_expmeter, :expmeter_rms => rms_expmeter, :expmeter_mean_winsor=>expmeter_mean_winsor, :expmeter_rms_winsor=>expmeter_rms_winsor, :expmeter_mean_red=>mean_expmeter_red, :expmeter_rms_red=>rms_expmeter_red, :expmeter_mean_green=>mean_expmeter_green, :expmeter_rms_green=>rms_expmeter_green, :expmeter_mean_blue=>mean_expmeter_blue, :expmeter_rms_blue=>rms_expmeter_blue )
    catch ex
       if typeof(f) == String
          @warn("*** Error extracting exposure meter data for " * string(f))
       else
          @warn("*** Error extracting exposure meter data for " * string(f.filename))
       end
       #return Dict(:expmeter_mean => missing, :expmeter_rms => missing, :expmeter_mean_winsor=>expmeter_mean_winsor, :expmeter_rms_winsor=>expmeter_rms_winsor, :expmeter_mean_red=>missing, :expmeter_mean_green=>missing, :expmeter_mean_blue=>missing )
       return Dict(:expmeter_mean => missing, :expmeter_rms => missing, :expmeter_mean_winsor=>missing, :expmeter_rms_winsor=>missing, :expmeter_mean_red=>missing, :expmeter_mean_green=>missing, :expmeter_mean_blue=>missing, :expmeter_rms_red=>missing, :expmeter_rms_green=>missing, :expmeter_rms_blue=>missing  )
    end
end

function add_metadata_from_fits!(df::DataFrame, fn::String)
    println("# Reading metadata from ", fn)
    metadata_keep = read_metadata(fn)
    #if metadata_keep[:target] != "Solar" return df end
    calibration_target_substrings = ["Etalon","LFC","ThArS","LDLS","FlatBB"]

   if any(map(ss->contains( metadata_keep[:target],ss),calibration_target_substrings))
      return df
   end


    for k in propertynames(df)
        if typeof(metadata_keep[k]) == Missing continue end
        if typeof(metadata_keep[k]) == Nothing
            metadata_keep[k] = missing
            continue
        end
        if (eltype(df[!,k]) <: Union{Real,Union{<:Real,Missing}}) && !(typeof(metadata_keep[k]) <: Union{Real,Union{<:Real,Missing}})
            metadata_keep[k] = missing
        end
    end
    #println("metadata_keep = ", metadata_keep)
    if haskey(metadata_keep,:airmass) && typeof(metadata_keep[:airmass]) == String
       metadata_keep[:airmass] = parse(Float64,metadata_keep[:airmass])
    end
    push!(df, metadata_keep, promote=true)
    return df
end


""" Read NEID data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(f::FITS, metadata::Dict{Symbol,Any}; normalization::Symbol = :raw, return_λ_obs::Bool=false)
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: Union{FITSIO.ImageHDU,FITSIO.TableHDU},2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["SCIWAVE"]]), FITSIO.read(f[img_idx["SCIFLUX"]]), FITSIO.read(f[img_idx["SCIVAR"]])

    if normalization == :raw
        metadata[:normalization] = :raw
    elseif normalization == :blaze
        blaze = read_blaze(f)
        flux ./= blaze
        var ./= blaze.^2
        metadata[:normalization] = :blaze
    else
        @error "# Reading data directly with normalization " * string(normalization) * " is not implemented."
    end

    if return_λ_obs
        spectrum = Spectra2DExtended(λ, copy(λ), flux, var, NEID2D(), metadata=metadata)
    else
        spectrum = Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
    end
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end


function read_data(f::FITS, metadata::Dict{Symbol,Any}, orders_to_read::AR; normalization::Symbol = :raw, return_λ_obs::Bool=false) where  { AR<:AbstractRange }
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: Union{FITSIO.ImageHDU,FITSIO.TableHDU},2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["SCIWAVE"]],:,orders_to_read), FITSIO.read(f[img_idx["SCIFLUX"]],:,orders_to_read), FITSIO.read(f[img_idx["SCIVAR"]],:,orders_to_read)

    if normalization == :raw
        metadata[:normalization] = :raw
    elseif normalization == :blaze
        blaze = read_blaze(f,orders_to_read)
        flux ./= blaze
        var ./= blaze.^2
        metadata[:normalization] = :blaze
    else
        @error "# Reading data directly with normalization " * string(normalization) * " is not implemented."
    end
    if return_λ_obs
        spectrum = Spectra2DExtended(λ, copy(λ), flux, var, NEID2D(), metadata=metadata)
    else
        spectrum = Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
    end
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end


function read_data(fn::String, metadata::Dict{Symbol,Any}; kwargs...)
    f = FITS(fn)
    read_data(f, metadata; kwargs...)
end

function read_data(fn::String; normalization::Symbol = :raw, kwargs...)
    f = FITS(fn)
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    read_data(f, metadata; normalization=normalization, kwargs...)
end

function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index}; kwargs...)
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn,metadata; kwargs...)
end

function read_data(fn::String, metadata::Dict{Symbol,Any}, orders_to_read::AR; kwargs...) where  { AR<:AbstractRange }
    f = FITS(fn)
    read_data(f, metadata, orders_to_read; kwargs...)
end

function read_data(fn::String, orders_to_read::AR; kwargs...) where  { AR<:AbstractRange }
    f = FITS(fn)
    metadata = read_metadata(f)
    read_data(f,metadata, orders_to_read; kwargs...)
end

function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index}, orders_to_read::AR; kwargs...) where  { AR<:AbstractRange }
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn, metadata, orders_to_read; kwargs...)
end

function apply_barycentric_correction!(λ::AbstractArray{T,1}, z::Real) where { T<:Real }
    # TODO: WRITE ME
end

""" Read blaze for NEID science fiber from FITS file."""
function read_blaze end

function read_blaze(fn::String)
    f = FITS(fn)
    read_blaze(f)
end

function read_blaze(fn::String, orders_to_read::AR ) where  { AR<:AbstractRange }
    f = FITS(fn)
    read_blaze(f,orders_to_read)
end

function read_blaze(f::FITS )
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: Union{FITSIO.ImageHDU,FITSIO.TableHDU},2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    blaze  = FITSIO.read(f[img_idx["SCIBLAZE"]])
    return blaze
end

function read_blaze(f::FITS, orders_to_read::AR ) where  { AR<:AbstractRange }
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: Union{FITSIO.ImageHDU,FITSIO.TableHDU},2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    blaze  = FITSIO.read(f[img_idx["SCIBLAZE"]],:,orders_to_read)
    return blaze
end


""" Read space delimited file with differential extinction corrections, interpolate to bjd's in df and insert into df[:,diff_ext_rv]. """
function read_differential_extinctions!(fn::String, df::DataFrame, df_time_col::Symbol = :bjd)
    df_diff_ext = CSV.read(fn, DataFrame, delim='\t')
    @assert any(isequal(:JD),propertynames(df_diff_ext))
    @assert any(isequal(:delta_vr),propertynames(df_diff_ext))
    diff_ext = LinearInterpolation(df_diff_ext.JD, df_diff_ext.delta_vr)
    df[!,:diff_ext_rv] = -diff_ext(df[!,df_time_col])/100 # cm/s -> m/s
    return df
end

""" Read normalization function for approximate solar continuum+blaze correction
TODO: Figure out mapping between these orders and what's stored in the FITS file before using. """
function read_solar_normalization(fn::String = joinpath(pkgdir(EchelleInstruments),"data","neid","NEID_normalization_function.csv") )
    CSV.read(fn, DataFrame)
end

""" Read calibrator data """
function read_cal_data(f::FITS, metadata::Dict{Symbol,Any}, orders_to_read::AR )  where { AR<:AbstractRange  }
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: Union{FITSIO.ImageHDU,FITSIO.TableHDU},2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["CALWAVE"]],:,orders_to_read), FITSIO.read(f[img_idx["CALFLUX"]],:,orders_to_read), FITSIO.read(f[img_idx["CALVAR"]],:,orders_to_read)
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end

function read_cal_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    read_cal_data(f, metadata)
end

function read_cal_data(fn::String)
    f = FITS(fn)
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    read_cal_data(f, metadata)
end
