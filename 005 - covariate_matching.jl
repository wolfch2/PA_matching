# https://discourse.julialang.org/t/elementwise-comparison-of-arrays-yield-a-x-element-bitarray-1-with-zeros-instead-of-bool/28700
# https://discourse.julialang.org/t/why-do-we-need-the-do-keyword/22236/10
# https://stackoverflow.com/questions/51459459/parallel-computing-in-julia-running-a-simple-for-loop-on-multiple-cores
# https://julialang.org/blog/2019/07/multithreading/
# https://stackoverflow.com/questions/29661315/vectorized-in-function-in-julia
# https://stackoverflow.com/questions/46411683/julia-check-if-elements-from-one-vector-are-within-another-vector
# https://discourse.julialang.org/t/how-can-i-rasterize-the-polygons-of-a-shapefile-in-julia/15089/4

using Distributed
addprocs(12)
@everywhere using ArchGDAL
using DataFrames
using Statistics
using Discretizers
using Shapefile
using FreqTables
using RCall
using Plots
using StatsBase
using CSV
using LibGEOS

@everywhere cd("/home/chrisgraywolf/shared/analysis/PA_matching/")

function full_dataset()
        keep = ArchGDAL.read("data_processed/rasters/cover.tif") do dataset
                band1 = ArchGDAL.getband(dataset,1)
                vec(ArchGDAL.read(band1) .>= 30) # only retain cells w/ at least 30% tree cover
        end
        ###
        rast_names = ["elev","slope",
                      "cover","loss","gain","cover_loss","lossyear",
                      "travel_time","pop_dens",
                      "countries","ecoregions","forest_biome",
                      "PAs","PAs_buffer",
                      "drivers",
                      "threatened","non_threatened"]
        all_data = pmap(WorkerPool(Array(2:4)), rast_names) do rast_name
                print(rast_name)
                rast_vec = ArchGDAL.read("data_processed/rasters/" * rast_name * ".tif") do dataset
                        band1 = ArchGDAL.getband(dataset,1)
                        vec(ArchGDAL.read(band1))[keep]
                end
        end
        full = DataFrame(all_data, Symbol.(rast_names))
        full[!,:pop_dens] = log.(1 .+ full[!,:pop_dens])
        full = full[(full.ecoregions .>= 0) .&
                    (full.countries .>= 0) .&
                    (full.drivers .>= 0) .&
                    (full.forest_biome .== 1) .&
                    (full.slope .>= 0) .&
                    (isnan.(full.cover_loss) .== false) .&
                    (isnan.(full.pop_dens) .== false),:]
        full = full[.!((full.PAs .< 0) .& (full.PAs_buffer .== 1)),:] # exclude points w/in ~10 km of PAs
        #full[:,:loss] = - (((full[:,:cover] .- full[:,:cover_loss]) ./ full[:,:cover]).^(1/18) .- 1)
        #select!(full, Not([:forest_biome, :cover_loss, :PAs_buffer]))
        return full
end

# calculate deforestation rate: overall, accounting for gain, and before/after if applicable (2002-2017)
function deforestation_calc(df, year, name_prefix="")
        C = sum(df[:,:cover])
        CL = sum(df[:,:cover_loss])
        CL_b = NaN
        CL_a = NaN
        Y_b = NaN
        Y_a = NaN
        L = sum(df[:,:loss])
        G = sum(df[:,:gain])
        if year in 2002:2017
                CL_b = sum(df[indexin(2000 .+ df.lossyear, 2001:year-1) .!= nothing,:cover_loss])
                CL_a = sum(df[indexin(2000 .+ df.lossyear, year+1:2018) .!= nothing,:cover_loss])
                Y_b = length(2001:year-1)
                Y_a = length(year+1:2018)
        end
        #
        out = DataFrame(loss = - (((C - CL)/C)^(1/18) - 1),
                        loss_before = - (((C - CL_b)/C)^(1/Y_b) - 1),
                        loss_after = - (((C - CL_a)/C)^(1/Y_a) - 1),
                        loss_gain = (L/18 - G/13)/C)
        rename!(out, [Symbol(name_prefix * String(x)) for x in names(out)])
        return out
end

# var contains names of columns to leave coarsen
function coarsen_df(df_original, vars=[], method="quantile", n_class=5)
        df = deepcopy(df_original)
        for col in vars # can prefix w/ Threads.@threads
                print(col)
                if method == "quantile" # DiscretizeUniformCount doesn't handle ties well (extremely slow and often fails...)
                        df[:,col] = encode(LinearDiscretizer(quantile(float(df[:,col]), Array(0:1/n_class:1))), float(df[:,col]))
                elseif method == "uniform"
                        df[:,col] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(n_class), float(df[:,col]))), float(df[:,col]))
                elseif method == "uniform_auto"
                        df[:,col] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(:sturges), float(df[:,col]))), float(df[:,col]))
                end
        end
        return df
end

function L1_imbalance(df)
        n_control = sum(df.treatment .== 0)
        n_treatment = sum(df.treatment .== 1)

        summary_gp = by(df, setdiff(names(df), ["treatment"])) do df_small # names(df) now returns list of strings, not symbols!
                abs(sum(df_small.treatment .== 0) / n_control - sum(df_small.treatment .== 1) / n_treatment)
        end

        return sum(summary_gp[:,:x1]) / 2
end

# perform 1-to-k coarsened exact matching (TODO - add response_var arg etc.).  That is, match,
# calculate loss rates for control pixels, then join to original treatment dataframe!
function CEM(data_treatment, data_control)
        ############# basic setup - add treatment variable and combine
        coarsen = [:elev, :slope, :cover, :travel_time, :pop_dens] # variables that must be coarsened        
        match = [:elev, :slope, :cover, :travel_time, :pop_dens, :countries, :ecoregions, :drivers] # variables to use for CEM
        # build combined dataset
        data_control[:,:treatment] .= 0
        data_treatment[:,:treatment] .= 1
        # https://stackoverflow.com/questions/51544317/breaking-change-on-vcat-when-columns-are-missing
        # https://discourse.julialang.org/t/vcat-does-not-append-dataframe-rows-as-it-did-in-0-10/9823/47
        for n in unique([names(data_control); names(data_treatment)]), df in [data_control,data_treatment]
                n in names(df) || (df[n] = NaN) # rbind.fill alternative - super ugly!
        end
        data_all = vcat(data_treatment, data_control)
        ############# aside: check L1 imbalance
        data_all_coarsened_auto = coarsen_df(data_all, coarsen, "uniform", 10)
        L1_imbalance(data_all_coarsened_auto[!,vcat(:treatment, match)]) # 0.879 (high imbalance)
        #############
        data_original = deepcopy(data_all[!,coarsen]) # needed for diagnostics; cover also needed for loss calc
        rename!(data_original, names(data_original) .* "_original")
        data_all = hcat(data_all, data_original)
        ############# now perform 1-k matching by iterating over the PA cells
        data_all_coarsened = coarsen_df(data_all, coarsen)
        data_control_coarsened = data_all_coarsened[data_all_coarsened.treatment .== 0,setdiff(names(data_all_coarsened),["PAs","STATUS_YR"])]
        data_treatment_coarsened = data_all_coarsened[data_all_coarsened.treatment .== 1,:]
        data_matched = join(data_treatment_coarsened[!,vcat(:PAs,:STATUS_YR,match)], data_control_coarsened, on = match, kind = :inner) # key matching step
        data_matched.cover = deepcopy(data_matched.cover_original) # non-coarsened version (for calculating loss rates)
        data_matched = by(data_matched, :PAs) do df_small # for each PA, want to calculate deforestation rate in the control pixels, etc.
                out = float.(DataFrame(df_small[1,:]))
                out = hcat(out, deforestation_calc(df_small, df_small.STATUS_YR[1], "Control_"))
                out[:,:n_control] .= length(df_small.loss)
                for var in names(data_original) # add control region means (for diagnostics)
                        out[:,var] = mean(df_small[:,var])
                end
                return out
        end # note: cover_original becomes spikey because we are averaging.
        # For example, see histogram(data_matched[data_matched.PAs .== 56624,:].cover_original)
        # cover must be in same 20% bracket, then we avg...
       ############# do a final join with original data_treatment and the control loss vars!
        out = join(data_treatment,  # has non-coarsened treatment covariates
                   data_matched[:,vcat([:PAs, :Control_loss, :Control_loss_before,
                                   :Control_loss_after, :Control_loss_gain, :n_control], Symbol.(names(data_original)))], # has matched info
                   on = :PAs,
                   kind = :left)
        return out
end

################################ build full dataset and remove recent PAs

full = full_dataset() # 26,254,884

PA_df = CSV.read("data_processed/PAs_tab.csv")
main_PAs = PA_df[(PA_df.MARINE .== 0) .& (PA_df.STATUS .!= "Proposed") .& (PA_df.GIS_AREA .>= 1),:]

# average covariates within PAs, adding year of establishment
data_treatment = full[indexin(full.PAs,main_PAs.ID) .!= nothing,:]
data_treatment[!,:STATUS_YR] = main_PAs.STATUS_YR[indexin(data_treatment.PAs,main_PAs.ID)]

# note: This gets cols. 2 and 3: out[1,[2,3]]
# but, this DOES NOT SET vals for these cols! out[1,[2,3]] = [1,2]
# instead, need to do broadcast: out[1,[2,3]] = [1,2]

data_treatment = by(data_treatment, :PAs) do df_small
        disc_vars = [:lossyear, :countries, :ecoregions, :PAs, :drivers]
        cont_vars = setdiff(names(df_small), disc_vars)
        #
        out = float.(DataFrame(df_small[1,:]))
        out[1,cont_vars] .= mean(convert(Array,df_small[:,cont_vars]), dims=1)[1,:]
        out[1,disc_vars] .= [mode(col) for col = eachcol(df_small[:,disc_vars]), row=1] #
        out = hcat(out, deforestation_calc(df_small, df_small.STATUS_YR[1], "PA_"))  
        return out
end

data_control = full[full.PAs .== minimum(full.PAs),:]

################################ run CEM and export

data_matched = CEM(data_treatment, data_control) # add control_loss

CSV.write("data_processed/data_matched.csv", data_matched)

