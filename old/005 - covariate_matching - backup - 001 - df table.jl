# https://discourse.julialang.org/t/elementwise-comparison-of-arrays-yield-a-x-element-bitarray-1-with-zeros-instead-of-bool/28700
# https://discourse.julialang.org/t/why-do-we-need-the-do-keyword/22236/10
# https://stackoverflow.com/questions/51459459/parallel-computing-in-julia-running-a-simple-for-loop-on-multiple-cores
# https://julialang.org/blog/2019/07/multithreading/
# https://stackoverflow.com/questions/29661315/vectorized-in-function-in-julia
# https://stackoverflow.com/questions/46411683/julia-check-if-elements-from-one-vector-are-within-another-vector
# https://discourse.julialang.org/t/how-can-i-rasterize-the-polygons-of-a-shapefile-in-julia/15089/4
#

using Distributed
addprocs(24)
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

@everywhere cd("/home/chrisgraywolf/analysis_desktop/PA_matching/")

function full_dataset()
        keep = ArchGDAL.read("data_processed/rasters/cover.tif") do dataset
                band1 = ArchGDAL.getband(dataset,1)
                vec(ArchGDAL.read(band1) .>= 30)
        end
        ###
        rast_names = ["elev","slope","cover","ecoregions","forest_biome","travel_time","PAs","PAs_buffer","cover_loss"]
        all_data = pmap(WorkerPool(Array(2:4)), rast_names) do rast_name
                print(rast_name)
                rast_vec = ArchGDAL.read("data_processed/rasters/" * rast_name * ".tif") do dataset
                        band1 = ArchGDAL.getband(dataset,1)
                        vec(ArchGDAL.read(band1))[keep]
                end
        end
        full = DataFrame(all_data, Symbol.(rast_names))
        full = full[(full.ecoregions .>= 0) .& (full.forest_biome .== 1) .& (full.slope .>= 0) .& (isnan.(full.cover_loss) .== false),:]
        full = full[.!((full.PAs .< 0) .& (full.PAs_buffer .== 1)),:] # exclude points w/in ~10 km of PAs
        full[:,:loss] = - (((full[:,:cover] .- full[:,:cover_loss]) ./ full[:,:cover]).^(1/18) .- 1)
        select!(full, Not([:forest_biome, :cover_loss, :PAs_buffer]))
        return full
end

# omit contains columns to leave alone (e.g., treatment variable, categorical vars.),
# manual is array containing names and breaks for specific columns (optional), method (either
# quantile or uniform) specifying procedure for dividing other variables into n_class
# groups. Example: manual = [("a",[1,3,5]), ("b",[2,4,6])]
function coarsen_df(df_original, omit=[], method="quantile", n_class=5, manual=[])
        df = deepcopy(df_original)
        for col in names(df) # can prefix w/ Threads.@threads
                print(col)
                if String(col) in omit
                        continue
                end
                if String(col) in [x[1] for x in manual]
                        breaks = manual[[x[1] == String(col) for x in manual]][1][2]
                        df[:,col] = encode(LinearDiscretizer(breaks), float(df[:,col]))
                elseif method == "quantile" # DiscretizeUniformCount doesn't handle ties well (extremely slow and often fails...)
                        # df[:,col] = encode(LinearDiscretizer(binedges(DiscretizeUniformCount(n_class), float(df[:,col]))), float(df[:,col]))
                        df[:,col] = encode(LinearDiscretizer(quantile(float(df[:,col]), Array(0:1/n_class:1))), float(df[:,col]))
                elseif method == "uniform"
                        df[:,col] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(n_class), float(df[:,col]))), float(df[:,col]))
                elseif method == "uniform_auto"
                        df[:,col] = encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(:sturges), float(df[:,col]))), float(df[:,col]))
                end
        end
        return df
end

# TODO: add a binary treatment_var argument...
function L1_imbalance(df, omit)
        n_control = sum(df.PAs .< 0)
        n_treatment = sum(df.PAs .>= 0)

        summary_gp = by(df, setdiff(names(df), Symbol.(omit))) do df_small
                abs(sum(df_small.PAs .< 0) / n_control - sum(df_small.PAs .>= 0) / n_treatment)
        end

        return sum(summary_gp[:,:x1]) / 2
end

# perform 1-to-k coarsened exact matching (TODO - add response_var arg etc.)
function CEM(data_treatment, data_control)
        omit = ["PAs","ecoregions","loss","treatment"]
        #
        data_control[:,:treatment] .= 0 # 20,686,942 rows
        data_treatment[:,:treatment] .= 1
        #
        data_all = vcat(data_treatment, data_control)
        data_all_coarsened_auto = coarsen_df(data_all, omit, "uniform", 10)
        L1_imbalance(data_all_coarsened_auto, omit) # 0.369
        #
        data_all_coarsened = coarsen_df(data_all, omit)
        #
        data_control_coarsened = data_all_coarsened[data_all_coarsened.treatment .== 0,:]
        data_control_coarsened = by(data_control_coarsened, setdiff(names(data_control_coarsened), [:loss])) do df_small
                out = DataFrame(df_small[1,:])
                out.loss .= mean(df_small.loss) # want average loss within groups
                out[:,:n_control] .= length(df_small.loss)
                return out
        end # 59,204 rows
        data_treatment_coarsened = data_all_coarsened[data_all_coarsened.treatment .== 1,:] # 16,256 rows
        # now join by covariates
        data_control_coarsened[:control_loss] = data_control_coarsened.loss
        select!(data_control_coarsened, Not([:PAs, :treatment, :loss]))
        data_joined = join(data_treatment_coarsened, data_control_coarsened, on = setdiff(names(data_control_coarsened), [:control_loss, :n_control]), kind = :inner)
        # plot(data_joined.loss, data_joined.control_loss, seriestype=:scatter)
        # cor(data_joined.loss, data_joined.control_loss)
        data_treatment = data_treatment[indexin(data_treatment.PAs,data_joined.PAs) .!= nothing,:]
        data_treatment = hcat(data_treatment, data_joined[indexin(data_treatment.PAs,data_joined.PAs),[:control_loss,:n_control]])
        # 
        return data_treatment # 12,859 rows
end

################################ build full dataset and remove recent PAs

full = full_dataset()

PA_table = Shapefile.Table("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp")
PA_df = DataFrame(PA_table)

# main PAs: Designated, area > 1km2, est. before 2000 // freqtable(PA_df.STATUS)
# PA_df.GIS_AREA[PA_df.NAME .== "Yellowstone National Park"] # units are km2
main_PAs = PA_df[(PA_df.MARINE .== "0") .& (PA_df.STATUS .!= "Proposed") .& (PA_df.GIS_AREA .>= 1) .& (PA_df.STATUS_YR .< 2000),:WDPAID]

# want to retain only main_PAs, not sure why .∈ is slow; compare with
# x = full.PAs
# y = main_PAs
# @rput x
# @rput y
# R"test = x %in% y"
# R"test = match(x,y)" # same thing..
# update: .∈ Ref(...) is awfully slow, but indexin is fine
# (probably sorts first and then does something smart that
# doesn't apply to the single element case)
# test = full.PAs .∈ Ref(intersect(full.PAs, main_PAs))

# average covariates within PAs
data_treatment = full[indexin(full.PAs,main_PAs) .!= nothing,:]
data_treatment = by(data_treatment, :PAs) do df_small
        out = float.(DataFrame(df_small[1,:]))
        cont = setdiff(names(df_small),[:ecoregions, :PAs]) # continuous variables to average over
        out[1,cont] = mean(convert(Array,df_small[:,cont]), dims=1)
        out.ecoregions .= mode(df_small.ecoregions)
        return out
end

data_control = full[full.PAs .== minimum(full.PAs),:]

################################ run CEM and export

data_matched = CEM(data_treatment, data_control) # add control_loss

CSV.write("data_processed/data_matched.csv", data_matched)

