##### Utility Functions for Focal Mechanism Stress Field Inversion Codes ########
##### Daniel Trugman, 2023
#####
##### Some codes adapted from the python package: https://github.com/ebeauce/ILSI 

##### Imports #####
using StatsBase
using LinearAlgebra
include("stress_utilities.jl")

### Function to Assemble Design Matrix for Inversion
# The normal vectors are in the coordinate system (x1, x2, x3):
#     x1: north
#     x2: west
#     x3: upward

function get_Gmat(norm1::Vector{Float32},norm2::Vector{Float32},norm3::Vector{Float32})
    
    # setup matrix
    nq = length(norm1)
    G = zeros(Float32, nq*3, 5)
    
    # loop over events, assemble matrix
    for i in eachindex(norm1)
        n1, n2, n3 = norm1[i], norm2[i], norm3[i]
        ii = 1 + 3*(i-1)
        G[ii + 0, 1] = n1 + n1 * n3^2 - n1^3
        G[ii + 0, 2] = n2 - 2.0f0 * n2 * n1^2
        G[ii + 0, 3] = n3 - 2.0f0 * n3 * n1^2
        G[ii + 0, 4] = n1 * n3^2 - n1 * n2^2
        G[ii + 0, 5] = -2.0f0 * n1 * n2 * n3
        G[ii + 1, 1] = n2 * n3^2 - n2 * n1^2
        G[ii + 1, 2] = n1 - 2.0f0 * n1 * n2^2
        G[ii + 1, 3] = -2.0f0 * n1 * n2 * n3
        G[ii + 1, 4] = n2 + n2 * n3^2 - n2^3
        G[ii + 1, 5] = n3 - 2.0f0 * n3 * n2^2
        G[ii + 2, 1] = n3^3 - n3 - n3 * n1^2
        G[ii + 2, 2] = -2.0f0 * n1 * n2 * n3
        G[ii + 2, 3] = n1 - 2.0f0 * n1 * n3^2
        G[ii + 2, 4] = n3^3 - n3 - n3 * n2^2
        G[ii + 2, 5] = n2 - 2.0f0 * n2 * n3^2
    end
    
    # return 
    return G
end

function get_Gmat(norm1::Vector{Float64},norm2::Vector{Float64},norm3::Vector{Float64})
    
    # setup matrix
    nq = length(norm1)
    G = zeros(Float64, nq*3, 5)
    
    # loop over events, assemble matrix
    for i in eachindex(norm1)
        n1, n2, n3 = norm1[i], norm2[i], norm3[i]
        ii = 1 + 3*(i-1)
        G[ii + 0, 1] = n1 + n1 * n3^2 - n1^3
        G[ii + 0, 2] = n2 - 2.0 * n2 * n1^2
        G[ii + 0, 3] = n3 - 2.0 * n3 * n1^2
        G[ii + 0, 4] = n1 * n3^2 - n1 * n2^2
        G[ii + 0, 5] = -2.0 * n1 * n2 * n3
        G[ii + 1, 1] = n2 * n3^2 - n2 * n1^2
        G[ii + 1, 2] = n1 - 2.0 * n1 * n2^2
        G[ii + 1, 3] = -2.0 * n1 * n2 * n3
        G[ii + 1, 4] = n2 + n2 * n3^2 - n2^3
        G[ii + 1, 5] = n3 - 2.0 * n3 * n2^2
        G[ii + 2, 1] = n3^3 - n3 - n3 * n1^2
        G[ii + 2, 2] = -2.0 * n1 * n2 * n3
        G[ii + 2, 3] = n1 - 2.0 * n1 * n3^2
        G[ii + 2, 4] = n3^3 - n3 - n3 * n2^2
        G[ii + 2, 5] = n2 - 2.0 * n2 * n3^2
    end
    
    # return 
    return G
end

### Function to Assemble Data Matrix for Inversion
# The slip vectors are in the coordinate system (x1, x2, x3):
#     x1: north
#     x2: west
#     x3: upward

function get_dvec(slip1::Vector{Float32},slip2::Vector{Float32},slip3::Vector{Float32})
    
    # setup vector
    nq = length(slip1)
    d = zeros(Float32,3*nq)
    
    # loop over events, assemble
    for i in eachindex(slip1)
        ii = 1 + 3*(i-1)
        d[ii] = slip1[i]
        d[ii+1] = slip2[i]
        d[ii+2] = slip3[i]
    end
    
    # return
    return d

end

function get_dvec(slip1::Vector{Float64},slip2::Vector{Float64},slip3::Vector{Float64})
    
    # setup vector
    nq = length(slip1)
    d = zeros(Float64,3*nq)
    
    # loop over events, assemble
    for i in eachindex(slip1)
        ii = 1 + 3*(i-1)
        d[ii] = slip1[i]
        d[ii+1] = slip2[i]
        d[ii+2] = slip3[i]
    end
    
    # return
    return d

end

### Function to Facilitate Linear Stress Inversion

# inputs: fault planes (strike, dip, rake)
# returns: stress tensor
function linear_stress_inversion(strike1::Vector{Float32},
    dip1::Vector{Float32},rake1::Vector{Float32})

    # compute normals and slips
    norm1, norm2, norm3 = get_normal(strike1,dip1)
    slip1, slip2, slip3 = get_slip(strike1, dip1, rake1)

    # assemble G-matrix and data vector
    G = get_Gmat(norm1, norm2, norm3)
    d = get_dvec(slip1, slip2, slip3)

    # invert, form stress tensor
    m = pinv(G) * d
    return assemble_stress_tensor(m)

end

# inputs: fault planes (strike, dip, rake)
# returns: stress tensor
function linear_stress_inversion(strike1::Vector{Float64},
    dip1::Vector{Float64},rake1::Vector{Float64})

    # compute normals and slips
    norm1, norm2, norm3 = get_normal(strike1,dip1)
    slip1, slip2, slip3 = get_slip(strike1, dip1, rake1)

    # assemble G-matrix and data vector
    G = get_Gmat(norm1, norm2, norm3)
    d = get_dvec(slip1, slip2, slip3)

    # invert, form stress tensor
    m = pinv(G) * d
    return assemble_stress_tensor(m)

end

### Random Fault Planes Version of the Above
# - Two versions, with and without aux planes pre-computed

# inputs: strike, dip, rake + number of trials
# returns: stress tensor
function linear_stress_inversion_randomfaults(strike1::Vector{Float32},
    dip1::Vector{Float32},rake1::Vector{Float32},ntrials::Int64)

    # initialize stress tensor
    S = zeros(Float32,3,3)

    # get aux planes
    strike2, dip2, rake2 = aux_plane(strike1,dip1,rake1)
    nq = length(strike1)

    # initialize
    strike, dip, rake = Vector{Float32}(undef,nq), Vector{Float32}(undef,nq), Vector{Float32}(undef,nq)

    # loop over trials
    for n in 1:ntrials
        ichoice = sample(1:2,nq) .== 1
        strike .= ifelse.(ichoice,strike1,strike2)
        dip .= ifelse.(ichoice,dip1,dip2)
        rake .= ifelse.(ichoice,rake1,rake2)
        S .+= linear_stress_inversion(strike,dip,rake)
    end

    # normalize and return
    #return S ./ sqrt(sum(S.^2))
    return S ./ norm(S)
    
end

# inputs: strike, dip, rake for both planes + number of trials
# returns: stress tensor
function linear_stress_inversion_randomfaults(strike1::Vector{Float32},
    dip1::Vector{Float32},rake1::Vector{Float32}, strike2::Vector{Float32},
    dip2::Vector{Float32},rake2::Vector{Float32}, ntrials::Int64)

    # initialize stress tensor
    S = zeros(Float32,3,3)

    # initialize
    nq = length(strike1)
    strike, dip, rake = Vector{Float32}(undef,nq), Vector{Float32}(undef,nq), Vector{Float32}(undef,nq)

    # loop over trials
    for n in 1:ntrials
        ichoice = sample(1:2,nq) .== 1
        strike .= ifelse.(ichoice,strike1,strike2)
        dip .= ifelse.(ichoice,dip1,dip2)
        rake .= ifelse.(ichoice,rake1,rake2)
        S .+= linear_stress_inversion(strike,dip,rake)
    end

    # normalize and return
    #return S ./ sqrt(sum(S.^2))
    return S ./ norm(S)
    
end

# inputs: strike, dip, rake + number of trials
# returns: stress tensor
function linear_stress_inversion_randomfaults(strike1::Vector{Float64},
    dip1::Vector{Float64},rake1::Vector{Float64},ntrials::Int64)

    # initialize stress tensor
    S = zeros(Float64,3,3)

    # get aux planes
    strike2, dip2, rake2 = aux_plane(strike1,dip1,rake1)
    nq = length(strike1)

    # initialize
    strike, dip, rake = Vector{Float64}(undef,nq), Vector{Float64}(undef,nq), Vector{Float64}(undef,nq)

    # loop over trials
    for n in 1:ntrials
        ichoice = sample(1:2,nq) .== 1
        strike .= ifelse.(ichoice,strike1,strike2)
        dip .= ifelse.(ichoice,dip1,dip2)
        rake .= ifelse.(ichoice,rake1,rake2)
        S .+= linear_stress_inversion(strike,dip,rake)
    end

    # normalize and return
    #return S ./ sqrt(sum(S.^2))
    return S ./ norm(S)
    
end

# inputs: strike, dip, rake for both planes + number of trials
# returns: stress tensor
function linear_stress_inversion_randomfaults(strike1::Vector{Float64},
    dip1::Vector{Float64},rake1::Vector{Float64}, strike2::Vector{Float64},
    dip2::Vector{Float64},rake2::Vector{Float64}, ntrials::Int64)

    # initialize stress tensor
    S = zeros(Float64,3,3)

    # initialize
    nq = length(strike1)
    strike, dip, rake = Vector{Float64}(undef,nq), Vector{Float64}(undef,nq), Vector{Float64}(undef,nq)

    # loop over trials
    for n in 1:ntrials
        ichoice = sample(1:2,nq) .== 1
        strike .= ifelse.(ichoice,strike1,strike2)
        dip .= ifelse.(ichoice,dip1,dip2)
        rake .= ifelse.(ichoice,rake1,rake2) 
        S .+= linear_stress_inversion(strike,dip,rake)
    end

    # normalize and return
    #return S ./ sqrt(sum(S.^2))
    return S ./ norm(S)
    
end

### Function for Iterative Stress Inversion ###

# inputs: strike, dip, rake, frictions, Niter, Ninit
# returns: 3x3 Stress tensor and optimized friction
function iterative_stress_inversion(
    strike1::Vector{Float32},dip1::Vector{Float32},rake1::Vector{Float32},
    strike2::Vector{Float32},dip2::Vector{Float32},rake2::Vector{Float32},
    frictions::Vector{Float32},Niter::Int64,Ninit::Int64)

    ## initialize stress tensor
    if Ninit == 0
        S0 = linear_stress_inversion(strike1,dip1,rake1)
    else
        S0 = linear_stress_inversion_randomfaults(strike1,dip1,rake1,strike2,dip2,rake2,Ninit)
    end
    S = zeros(Float32,3,3)

    ## store results for each friction
    nfric = length(frictions)
    istabs = Vector{Float32}(undef,nfric)


    ## loop over frictions
    for (ii, friction) in enumerate(frictions)
        
        # initialize stress tensor for this friction
        S.=S0
        
        # loop over iterations
        for jj in 1:Niter
            
            # get principal stresses
            spvecs, spvals, R = principal_stresses(S)
            
            # optimize fault planes, given S
            I, strike, dip, rake = compute_instability(
                spvecs,R,friction,strike1,dip1,rake1,strike2, dip2,rake2)
            
            # compute S, given fault planes
            S .= linear_stress_inversion(strike,dip,rake)
            
            # compute final instability 
            if jj == Niter
                istabs[ii] = mean(I)
            end
            
        end
        
    end

    ## determine best friction
    ifric = argmax(istabs)
    friction = frictions[ifric]

    ## final run with optimized friction
    S .= S0
    for jj in 1:Niter

        # get principal stresses
        spvecs, spvals, R = principal_stresses(S)

        # optimize fault planes, given S
        I, strike, dip, rake = compute_instability(
            spvecs,R,friction,strike1,dip1,rake1,strike2, dip2,rake2)

        # compute S, given fault planes
        S .= linear_stress_inversion(strike,dip,rake)

    end

    ## return
    return S, friction

end

# inputs: strike, dip, rake, frictions, Niter, Ninit
# returns: 3x3 Stress tensor and optimized friction
function iterative_stress_inversion(
        strike1::Vector{Float64},dip1::Vector{Float64},rake1::Vector{Float64},
        strike2::Vector{Float64},dip2::Vector{Float64},rake2::Vector{Float64},
        frictions::Vector{Float64},Niter::Int64,Ninit::Int64)
    
    ## initialize stress tensor
    if Ninit == 0
        S0 = linear_stress_inversion(strike1,dip1,rake1)
    else
        S0 = linear_stress_inversion_randomfaults(strike1,dip1,rake1,strike2,dip2,rake2,Ninit)
    end
    S = zeros(Float64,3,3)
    
    ## store results for each friction
    nfric = length(frictions)
    istabs = Vector{Float64}(undef,nfric)
    
    
    ## loop over frictions
    for (ii, friction) in enumerate(frictions)
        
        # initialize stress tensor for this friction
        S.=S0
        
        # loop over iterations
        for jj in 1:Niter
            
            # get principal stresses
            spvecs, spvals, R = principal_stresses(S)
            
            # optimize fault planes, given S
            I, strike, dip, rake = compute_instability(
                spvecs,R,friction,strike1,dip1,rake1,strike2, dip2,rake2)
            
            # compute S, given fault planes
            S .= linear_stress_inversion(strike,dip,rake)
            
            # compute final instability 
            if jj == Niter
                istabs[ii] = mean(I)
            end
            
        end
        
    end
    
    ## determine best friction
    ifric = argmax(istabs)
    friction = frictions[ifric]

    ## final run with optimized friction
    S .= S0
    for jj in 1:Niter

        # get principal stresses
        spvecs, spvals, R = principal_stresses(S)

        # optimize fault planes, given S
        I, strike, dip, rake = compute_instability(
            spvecs,R,friction,strike1,dip1,rake1,strike2, dip2,rake2)

        # compute S, given fault planes
        S .= linear_stress_inversion(strike,dip,rake)

    end
    
    ## return
    return S, friction
    
end

#### Similar to the Above but with Mechanism Matrix (N x 6: s1, d1, r1, s2, d2, r2)

function iterative_stress_inversion(
    mechMAT::Matrix{Float32},
    frictions::Vector{Float32},Niter::Int64,Ninit::Int64)

    ## initialize stress tensor
    if Ninit == 0
        S0 = linear_stress_inversion(mechMAT[:,1],mechMAT[:,2],mechMAT[:,3])
    else
        S0 = linear_stress_inversion_randomfaults(
            mechMAT[:,1],mechMAT[:,2],mechMAT[:,3],
            mechMAT[:,4],mechMAT[:,5],mechMAT[:,6],Ninit)
    end
    S = zeros(Float32,3,3)

    ## store results for each friction
    nfric = length(frictions)
    istabs = Vector{Float32}(undef,nfric)


    ## loop over frictions
    for (ii, friction) in enumerate(frictions)
        
        # initialize stress tensor for this friction
        S.=S0
        
        # loop over iterations
        for jj in 1:Niter
            
            # get principal stresses
            spvecs, spvals, R = principal_stresses(S)
            
            # optimize fault planes, given S
            I, strike, dip, rake = compute_instability(
                spvecs,R,friction,mechMAT[:,1],mechMAT[:,2],mechMAT[:,3],
                mechMAT[:,4],mechMAT[:,5],mechMAT[:,6])
            
            # compute S, given fault planes
            S .= linear_stress_inversion(mechMAT[:,1],mechMAT[:,2],mechMAT[:,3])
            
            # compute final instability 
            if jj == Niter
                istabs[ii] = mean(I)
            end
            
        end
        
    end

    ## determine best friction
    ifric = argmax(istabs)
    friction = frictions[ifric]

    ## final run with optimized friction
    S .= S0
    for jj in 1:Niter

        # get principal stresses
        spvecs, spvals, R = principal_stresses(S)

        # optimize fault planes, given S
        I, strike, dip, rake = compute_instability(
            spvecs,R,friction,mechMAT[:,1],mechMAT[:,2],mechMAT[:,3],
            mechMAT[:,4],mechMAT[:,5],mechMAT[:,6])

        # compute S, given fault planes
        S .= linear_stress_inversion(mechMAT[:,1],mechMAT[:,2],mechMAT[:,3])

    end

    ## return
    return S, friction

end

function iterative_stress_inversion(
    mechMAT::Matrix{Float64},
    frictions::Vector{Float64},Niter::Int64,Ninit::Int64)

    ## initialize stress tensor
    if Ninit == 0
        S0 = linear_stress_inversion(mechMAT[:,1],mechMAT[:,2],mechMAT[:,3])
    else
        S0 = linear_stress_inversion_randomfaults(
            mechMAT[:,1],mechMAT[:,2],mechMAT[:,3],
            mechMAT[:,4],mechMAT[:,5],mechMAT[:,6],Ninit)
    end
    S = zeros(Float64,3,3)

    ## store results for each friction
    nfric = length(frictions)
    istabs = Vector{Float64}(undef,nfric)


    ## loop over frictions
    for (ii, friction) in enumerate(frictions)
        
        # initialize stress tensor for this friction
        S.=S0
        
        # loop over iterations
        for jj in 1:Niter
            
            # get principal stresses
            spvecs, spvals, R = principal_stresses(S)
            
            # optimize fault planes, given S
            I, strike, dip, rake = compute_instability(
                spvecs,R,friction,mechMAT[:,1],mechMAT[:,2],mechMAT[:,3],
                mechMAT[:,4],mechMAT[:,5],mechMAT[:,6])
            
            # compute S, given fault planes
            S .= linear_stress_inversion(mechMAT[:,1],mechMAT[:,2],mechMAT[:,3])
            
            # compute final instability 
            if jj == Niter
                istabs[ii] = mean(I)
            end
            
        end
        
    end

    ## determine best friction
    ifric = argmax(istabs)
    friction = frictions[ifric]

    ## final run with optimized friction
    S .= S0
    for jj in 1:Niter

        # get principal stresses
        spvecs, spvals, R = principal_stresses(S)

        # optimize fault planes, given S
        I, strike, dip, rake = compute_instability(
            spvecs,R,friction,mechMAT[:,1],mechMAT[:,2],mechMAT[:,3],
            mechMAT[:,4],mechMAT[:,5],mechMAT[:,6])

        # compute S, given fault planes
        S .= linear_stress_inversion(mechMAT[:,1],mechMAT[:,2],mechMAT[:,3])

    end

    ## return
    return S, friction

end


### Bootstrapping interative inversion: returns two quantile values only: Float32
#   - Inputs are same as iterative stress inversion, plus Nboot=Num bootstrap resamples, Qpts: two quantiles, 
#     and Sazi0: azimuths from non-bootstrap version
function bootstrap_iterative_inversion(strike1::Vector{Float32},dip1::Vector{Float32},rake1::Vector{Float32},
    strike2::Vector{Float32},dip2::Vector{Float32},rake2::Vector{Float32},
    frictions::Vector{Float32},Niter::Int64,Ninit::Int64,Nboot::Int64,Qpts::Vector{Float32},Sazi0::Vector{Float32})

    # allocated memory for bootstrap distribution of parameters
    fricB, misfitB, shapeB = Vector{Float32}(undef,Nboot), Vector{Float32}(undef,Nboot), Vector{Float32}(undef,Nboot)
    SaziB, SpluB = Matrix{Float32}(undef,Nboot,3), Matrix{Float32}(undef,Nboot,3)

    # sampling vector
    Ndat = length(strike1)
    idx = Vector{Int32}(1:Ndat)

    # loop over bootstraps resamplings
    Threads.@threads for ib in 1:Nboot

        # get sampling indices
        isamp = sample(idx,Ndat,replace=true)

        # run inversion
        S, friction = iterative_stress_inversion(
            strike1[isamp],dip1[isamp],rake1[isamp],
            strike2[isamp],dip2[isamp],rake2[isamp],
            frictions,Niter,Ninit)
        fricB[ib] = friction

        # compute outputs
        # eigen decomposition for principal stresses
        spvecs, _, R = principal_stresses(S)
        shapeB[ib] = R

        # optimize fault planes, given S
        _, strike, dip, rake = compute_instability(
            spvecs,R,friction,
            strike1[isamp],dip1[isamp],rake1[isamp],
            strike2[isamp],dip2[isamp],rake2[isamp]) 

        # get azimuths and plunges
        for nn in 1:3

            # solve for angles
            SaziB[ib,nn], SpluB[ib,nn] = get_azi_plunge(spvecs[:,nn])

            # special section to prevent wraparound weirdness
            if SaziB[ib,nn]-Sazi0[nn]>=270.0f0 
                SaziB[ib,nn]-=360.0f0
            elseif SaziB[ib,nn]-Sazi0[nn]>=90.0f0 
                SaziB[ib,nn]-=180.0f0
            elseif SaziB[ib,nn]-Sazi0[nn]<=-270.0f0
                SaziB[ib,nn]+=360.0f0
            elseif SaziB[ib,nn]-Sazi0[nn]<=-90.0f0
                SaziB[ib,nn]+=180.0f0
            end
        end

        # compute misfit and update count
        misfitB[ib] = mean_angular_misfit(S,strike,dip,rake)

    end

    # extract quantiles
    fricQ, shapeQ, misfitQ = quantile(fricB,Qpts), quantile(shapeB,Qpts), quantile(misfitB,Qpts)
    s1aQ, s2aQ, s3aQ = quantile(SaziB[:,1],Qpts), quantile(SaziB[:,2],Qpts), quantile(SaziB[:,3],Qpts) # differential
    s1pQ, s2pQ, s3pQ = quantile(SpluB[:,1],Qpts), quantile(SpluB[:,2],Qpts), quantile(SpluB[:,3],Qpts)

    # return from this
    return fricQ, shapeQ, misfitQ, [s1aQ s2aQ s3aQ], [s1pQ s2pQ s3pQ]

end

### Bootstrapping interative inversion: returns two quantile values only: Float64
#   - Inputs are same as iterative stress inversion, plus Nboot=Num bootstrap resamples, Qpts: two quantiles, 
#     and Sazi0: azimuths from non-bootstrap version
function bootstrap_iterative_inversion(strike1::Vector{Float64},dip1::Vector{Float64},rake1::Vector{Float64},
    strike2::Vector{Float64},dip2::Vector{Float64},rake2::Vector{Float64},
    frictions::Vector{Float64},Niter::Int64,Ninit::Int64,Nboot::Int64,Qpts::Vector{Float32},Sazi0::Vector{Float64})

    # allocated memory for bootstrap distribution of parameters
    fricB, misfitB, shapeB = Vector{Float64}(undef,Nboot), Vector{Float64}(undef,Nboot), Vector{Float64}(undef,Nboot)
    SaziB, SpluB = Matrix{Float64}(undef,Nboot,3), Matrix{Float64}(undef,Nboot,3)

    # sampling vector
    Ndat = length(strike1)
    idx = Vector{Int32}(1:Ndat)

    # loop over bootstraps resamplings
    Threads.@threads for ib in 1:Nboot

        # get sampling indices
        isamp = sample(idx,Ndat,replace=true)

        # run inversion
        S, friction = iterative_stress_inversion(
            strike1[isamp],dip1[isamp],rake1[isamp],
            strike2[isamp],dip2[isamp],rake2[isamp],
            frictions,Niter,Ninit)
        fricB[ib] = friction

        # compute outputs
        # eigen decomposition for principal stresses
        spvecs, _, R = principal_stresses(S)
        shapeB[ib] = R

        # optimize fault planes, given S
        _, strike, dip, rake = compute_instability(
            spvecs,R,friction,
            strike1[isamp],dip1[isamp],rake1[isamp],
            strike2[isamp],dip2[isamp],rake2[isamp]) 

        # get azimuths and plunges
        for nn in 1:3

            # compute initial angles
            SaziB[ib,nn], SpluB[ib,nn] = get_azi_plunge(spvecs[:,nn])

            # special section to prevent wraparound weirdness
            if SaziB[ib,nn]-Sazi0[nn]>=270.0 
                SaziB[ib,nn]-=360.0
            elseif SaziB[ib,nn]-Sazi0[nn]>=90.0 
                SaziB[ib,nn]-=180.0
            elseif SaziB[ib,nn]-Sazi0[nn]<=-270.0
                SaziB[ib,nn]+=360.0
            elseif SaziB[ib,nn]-Sazi0[nn]<=-90.0
                SaziB[ib,nn]+=180.0
            end
        end

        # compute misfit and update count
        misfitB[ib] = mean_angular_misfit(S,strike,dip,rake)

    end

    # extract quantiles
    fricQ, shapeQ, misfitQ = quantile(fricB,Qpts), quantile(shapeB,Qpts), quantile(misfitB,Qpts)
    s1aQ, s2aQ, s3aQ = quantile(SaziB[:,1],Qpts), quantile(SaziB[:,2],Qpts), quantile(SaziB[:,3],Qpts)
    s1pQ, s2pQ, s3pQ = quantile(SpluB[:,1],Qpts), quantile(SpluB[:,2],Qpts), quantile(SpluB[:,3],Qpts)

    # return from this
    return fricQ, shapeQ, misfitQ, [s1aQ s2aQ s3aQ], [s1pQ s2pQ s3pQ]

end