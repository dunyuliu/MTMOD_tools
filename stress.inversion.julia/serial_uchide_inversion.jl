######## Julia Script to Test Space-Time Inversion of Uchide et al. 2020/22 Mechanisms ##########
###      Daniel Trugman, 2023

######## Packages ########

# load Packages
using StatsBase
using LinearAlgebra
using Printf
using DelimitedFiles
using DataFrames
using Dates
using LazyGrids
using CSV

# and stress Functions
include("stress_inversion.jl")
include("stress_utilities.jl")

######### Define Program Parameters ##########

# data
# this data_dir is for mtmod docker image. 
# need to change it for other systems.
const data_dir = "/home/MTMOD_tools/stress.inversion.julia/" 
const mech_file = data_dir * "meca_Japan_Uchide2022.txt"

# floating point precision for calculations
const fptype = Float32
if fptype == Float64
    fpstr = "f64"
else
    fpstr = "f32"
end

# binning in X and T (bins may overlap if stride != bin)
const binX = binY = 0.5 # degrees
const strideX = strideY = 0.5*binX
const binT = 0.5 # years
const strideT = 0.5*binT 

# space time bounds
const lon0 = 128.0
const lon1 = 145.5
const lat0 = 27.5
const lat1 = 45.5
const year0 = 2002.0
const year1 = 2021.0

# parameters for inversion
const minsum = 200 # minimum number of points in a spatial bin, skip otherwise
const mincalc = 20 # minimum number of points in space-time bin to run inversion
const fmin = 0.2 # minimum friction
const fmax = 0.8 # maximum friction
const dfric = 0.05 # spacing between friction
const Niter = 5 # number of iterations to select best planes
const Ninit = 5 # number of initializations for stress tensor, does not do the bootstrapping for Ninit = 0

# outfile
out_dir = data_dir * "out1/"
if !(isdir(out_dir)); mkpath(out_dir); end
outfile = @sprintf("stress_inversion_binX-%.2f_strideX-%.2f_binY-%.2f_strideY-%.2f_binT-%.2f_strideT-%.2f.%s.txt",
    binX,strideX,binY,strideY,binT,strideT,fpstr)


###### Define Local Functions #########

# function to load Uchide et al. dataset into a DataFrame
function load_uchide_dataframe(mech_file, fptype)

    # define columns
    cols = ["lon","lat","depth","strike1","dip1","rake1",
    "mag","x","y","qual","nsta","fptype","otime"]

    # load dataframe
    df = DataFrame(readdlm(mech_file,Any,skipstart=2),cols)

    # subset
    cols = ["lon","lat","strike1","dip1","rake1","otime"]
    select!(df,cols)

    # datatype conversion
    for col in cols[1:end-1]
        df[!,col] = convert.(fptype,df[!,col])
    end

    # compute aux plane
    s2, d2, r2 = aux_plane(df.strike1,df.dip1,df.rake1)
    df[!,"strike2"] .= s2
    df[!,"dip2"] .= d2
    df[!,"rake2"] .= r2

    # reorganize new columns
    cols = ["lon","lat","strike1","dip1","rake1","strike2","dip2","rake2","otime"]
    select!(df,cols)

    # datetime calculations: Japan is 9 hours ahead of UTC
    dfmt = DateFormat("yyyy-mm-ddTHH:MM:SS:ssssss")
    df[!,"otime"] = [DateTime(x,dfmt) - Hour(9) for x in df.otime]

    # now calculate decimal year
    yptype = Float64
    years=Vector{yptype}(undef,nrow(df))
    millisec_per_day = 24*60*60*1000
    for (ii,tt) in enumerate(df.otime)
        nsec = convert(yptype,millisec_per_day*daysinyear(tt))
        year = convert(yptype,Year(tt).value)
        dsec = convert(yptype,(tt - DateTime(Year(tt))).value)
        years[ii] = year + dsec/nsec
    end
    df[!,"year"] = convert.(fptype,years)

    # return on necessary columns
    select!(df,["lon","lat","strike1","dip1","rake1","strike2","dip2","rake2","year"])
    return df
end

#################################################

################ Main Program ###################

#### Setup for Inversion ####

# load dataframe
df = load_uchide_dataframe(mech_file,fptype)
println(first(df,10))
println()
println(last(df,10))
println()

# extract arrays,
qlon, qlat, qtime = df.lon, df.lat, df.year
strike1, dip1, rake1 = df.strike1, df.dip1, df.rake1
strike2, dip2, rake2 = df.strike2, df.dip2, df.rake2

# define friction vector
const frictions = Vector{fptype}(fmin:dfric:fmax)

# define space-time bins
const lonC = Vector{fptype}(lon0:strideX:lon1)
println("Longitudes:")
println(lonC,"\n")
const latC = Vector{fptype}(lat0:strideY:lat1)
println("Latitudes:")
println(latC,"\n")
const timeC = Vector{fptype}(year0:strideT:year1)
println("Times:")
println(timeC,"\n")

# allocate memory for results
nI, nJ, nK = length(lonC), length(latC), length(timeC)
nTOT = nI*nJ*nK
@printf("Total bins: %d x %d x %d = %d\n",nI,nJ,nK,nI*nJ*nK)
if fptype == Float32
    RatMat = fill(NaN32,nI,nJ,nK)
    AziMat = fill(NaN32,nI,nJ,nK,3)
    PluMat = fill(NaN32,nI,nJ,nK,3)
    FriMat = fill(NaN32,nI,nJ,nK)
    MisMat = fill(NaN32,nI,nJ,nK)
else
    RatMat = fill(NaN64,nI,nJ,nK)
    AziMat = fill(NaN64,nI,nJ,nK,3)
    PluMat = fill(NaN64,nI,nJ,nK,3)
    FriMat = fill(NaN64,nI,nJ,nK)
    MisMat = fill(NaN64,nI,nJ,nK)
end
NumMat = zeros(Int64,nI,nJ,nK)

##### Run Inversions ####

# loop over longitude bins
@time for ii in 1:nI
    
    # select data
    lonX = lonC[ii]
    ix = (abs.(qlon .- lonX) .<= binX)
    
    # loop over latitude bins
    @printf("Working on Longitude bin: %d/%d...\n",ii,nI)
    @time for (jj,latX) in enumerate(latC)
        
        # select data
        ixy = (ix) .& (abs.(qlat .- latX) .<= binY)
        nxy = sum(ixy)
        @printf("Latitude: %d/%d, ndata=%d",jj,nJ,nxy)
        
        # skip if needed
        if nxy < minsum;
            @printf(" (skipping)\n") 
            continue
        else
            @printf("\n")
        end
        
        # loop over temporal bins
        for (tt,tC) in enumerate(timeC)
            
            # select data
            ixyt = (ixy) .& (abs.(qtime .- tC) .<= binT)
            count = sum(ixyt)
            if sum(ixyt) < mincalc; continue; end

            # run inversion
            S, friction = iterative_stress_inversion(
                strike1[ixyt],dip1[ixyt],rake1[ixyt],
                strike2[ixyt],dip2[ixyt],rake2[ixyt],
                frictions,Niter,Ninit)
            FriMat[ii,jj,tt] = friction

            # compute outputs
            # eigen decomposition for principal stresses
            spvecs, _, R = principal_stresses(S)
            RatMat[ii,jj,tt] = R

            # optimize fault planes, given S
            _, strike, dip, rake = compute_instability(
               spvecs,R,friction,
                strike1[ixyt],dip1[ixyt],rake1[ixyt],
                strike2[ixyt],dip2[ixyt],rake2[ixyt]) 

            # get azimuths and plunges
            for nn in 1:3
                AziMat[ii,jj,tt,nn], PluMat[ii,jj,tt,nn] = get_azi_plunge(spvecs[:,nn])
            end

            # compute misfit and update count
            MisMat[ii,jj,tt] = mean_angular_misfit(S,strike,dip,rake)
            NumMat[ii,jj,tt] = count
        end
    end
    @printf("Done with Longitude bin: %d/%d...\n\n\n",ii,nI)
    #if ii ==6; break; end # good stopping point for tests
end

### Format Output results
println("\n\nINVERSIONS COMPLETE")
iout = NumMat .> 0 # counts bins processed
@printf("Bins processed: %d/%d\n",sum(iout),nTOT)
println("Formatting output...")

# define lon, lat, time grids
xG, yG, tG = ndgrid_array(lonC,latC,timeC) 

# unpack principal stresses into columnizable vectors
s1a, s2a, s3a = AziMat[:,:,:,1], AziMat[:,:,:,2], AziMat[:,:,:,3]
s1p, s2p, s3p = PluMat[:,:,:,1], PluMat[:,:,:,2], PluMat[:,:,:,3]

# create dataframe
odf = DataFrame(    
    lon=vec(xG[iout]),
    lat=vec(yG[iout]),
    time=vec(tG[iout]),
    count=vec(NumMat[iout]),
    friction=vec(FriMat[iout]),
    misfit=vec(MisMat[iout]),
    shape=vec(RatMat[iout]),
    s1a=vec(s1a[iout]),
    s1p=vec(s1p[iout]),
    s2a=vec(s2a[iout]),
    s2p=vec(s2p[iout]),
    s3a=vec(s3a[iout]),
    s3p=vec(s3p[iout]),
)
println(first(odf,20))
println()

# output file
println("Output file: ",outfile)
CSV.write(out_dir*outfile, odf, delim=" ")
