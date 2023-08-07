##### Utility Functions for Focal Mechanism Stress Field Inversion Codes ########
##### Daniel Trugman, 2023
#####
##### Some codes adapted from the python package: https://github.com/ebeauce/ILSI 

#### Imports ####
using StatsBase
using LinearAlgebra

### Functions to Compute Normal and Slip Vectors
#
# The vectors are in the coordinate system (x1, x2, x3):
#     x1: north
#     x2: west
#     x3: upward

# Returns inward normal for hanging wall, following Wysession and Stein
function get_normal(strike::Vector{Float32},dip::Vector{Float32})
    norm1 = -sind.(dip) .* sind.(strike)
    norm2 = -sind.(dip) .* cosd.(strike)
    norm3 = cosd.(dip)
    return norm1, norm2, norm3
end

# Returns inward normal for hanging wall, following Wysession and Stein
function get_normal(strike::Vector{Float64},dip::Vector{Float64})
    norm1 = -sind.(dip) .* sind.(strike)
    norm2 = -sind.(dip) .* cosd.(strike)
    norm3 = cosd.(dip)
    return norm1, norm2, norm3
end

# Returns slip vector of the hanging wall w.r.t the footwall
function get_slip(strike::Vector{Float32},dip::Vector{Float32},rake::Vector{Float32})
    slip1 = cosd.(rake) .* cosd.(strike) .+ sind.(rake) .* cosd.(dip) .* sind.(strike)
    slip2 = -cosd.(rake) .* sind.(strike) .+ sind.(rake) .* cosd.(dip) .* cosd.(strike)
    slip3 = sind.(rake) .* sind.(dip)
    return slip1, slip2, slip3
end


# Returns slip vector of the hanging wall w.r.t the footwall
function get_slip(strike::Vector{Float64},dip::Vector{Float64},rake::Vector{Float64})
    slip1 = cosd.(rake) .* cosd.(strike) .+ sind.(rake) .* cosd.(dip) .* sind.(strike)
    slip2 = -cosd.(rake) .* sind.(strike) .+ sind.(rake) .* cosd.(dip) .* cosd.(strike)
    slip3 = sind.(rake) .* sind.(dip)
    return slip1, slip2, slip3
end


### Functions for Aux Plane Calculations
#
#    Taken from ObsPy: https://docs.obspy.org/_modules/obspy/imaging/beachball.html
#     Originally Adapted from MATLAB script
#     `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
#     written by Andy Michael, Chen Ji and Oliver Boyd.

## Auxiliary Plane Given Strike, Dip, Rake 
function aux_plane(s0::Vector{Float32}, d1::Vector{Float32}, r1::Vector{Float32})
    
    # redefine strike direction
    s1 = s0 .+ 90.0f0
    
    # slick vector in plane 1
    sl1 = -cosd.(r1) .* cosd.(s1) .- sind.(r1) .* sind.(s1) .* cosd.(d1)
    sl2 = cosd.(r1) .* sind.(s1) .- sind.(r1) .* cosd.(s1) .* cosd.(d1)
    sl3 = sind.(r1) .* sind.(d1)
    
    # get strike and dip
    strike, dip = strike_dip(sl2,sl1,sl3)

    # normal vector to plane 1
    n1 = sind.(s1) .* sind.(d1)  
    n2 = cosd.(s1) .* sind.(d1)
    # n3 = cosd(s2)
    
    # strike vector of plane 2
    h1 = -sl2  
    h2 = sl1
    # note h3=0 always so we leave it out
    
    # recompute z
    z = h1 .* n1 .+ h2 .* n2
    z ./= sqrt.(h1 .* h1 .+ h2 .* h2)
    
    # correct for roundoff error (both methods ok?)
    z = ifelse.(abs.(z).>=1.0,sign.(z),z)
    #z = ifelse.(isapprox.(abs.(z),1.0),sign.(z),z)
    
    # calculate rake
    z = acosd.(z)
    rake = ifelse.(sl3 .> 0.0f0, z, -z)
    
    # return
    return strike, dip, rake

end

## Auxiliary Plane Given Strike, Dip, Rake 
function aux_plane(s0::Vector{Float64}, d1::Vector{Float64}, r1::Vector{Float64})
    
    # redefine strike direction
    s1 = s0 .+ 90.0
    
    # slick vector in plane 1
    sl1 = -cosd.(r1) .* cosd.(s1) .- sind.(r1) .* sind.(s1) .* cosd.(d1)
    sl2 = cosd.(r1) .* sind.(s1) .- sind.(r1) .* cosd.(s1) .* cosd.(d1)
    sl3 = sind.(r1) .* sind.(d1)
    
    # get strike and dip
    strike, dip = strike_dip(sl2,sl1,sl3)

    # normal vector to plane 1
    n1 = sind.(s1) .* sind.(d1)  
    n2 = cosd.(s1) .* sind.(d1)
    # n3 = cosd(s2)
    
    # strike vector of plane 2
    h1 = -sl2  
    h2 = sl1
    # note h3=0 always so we leave it out
    
    # recompute z and normalize
    z = h1 .* n1 .+ h2 .* n2
    z ./= sqrt.(h1.^2 .+ h2.^2)
    
    # correct for roundoff error (both methods ok?)
    z = ifelse.(abs.(z).>=1.0,sign.(z),z)
    #z = ifelse.(isapprox.(abs.(z),1.0),sign.(z),z)
    
    # calculate rake
    z = acosd.(z)
    rake = ifelse.(sl3 .> 0.0, z, -z)
    
    # return
    return strike, dip, rake

end
    
## Strike and Dip Angles 
function strike_dip(n::Vector{Float32},e::Vector{Float32},u::Vector{Float32})
    
    # flip negatives
    iflip = u .< 0.0f0
    n = ifelse.(iflip,-n,n)
    e = ifelse.(iflip,-e,e)
    u = ifelse.(iflip,-u,u)
    
    # compute strike
    strike = atand.(e,n) .- 90.0f0
    strike = ifelse.(strike .< 0.0f0, strike.+360.0f0, strike)
    
    # compute dip
    x = sqrt.(n.^2 .+ e.^2)
    dip = atand.(x,u)
    
    # return
    return strike, dip
    
end

## Strike and Dip Angles 
function strike_dip(n::Vector{Float64},e::Vector{Float64},u::Vector{Float64})
    
    # flip negatives
    iflip = u .< 0.0
    n = ifelse.(iflip,-n,n)
    e = ifelse.(iflip,-e,e)
    u = ifelse.(iflip,-u,u)
    
    # compute strike
    strike = atand.(e,n) .- 90.0
    strike = ifelse.(strike .< 0.0, strike.+360.0, strike)
    
    # compute dip
    x = sqrt.(n.^2 .+ e.^2)
    dip = atand.(x,u)
    
    # return
    return strike, dip
    
end

### Given Stress Vector from Inversion, Assemble 3x3 Symmetric Tensor
#   Note, we are following Beauce's codes so s = [S11,S12,S13,S22,S23]
function assemble_stress_tensor(svec::Vector{Float32})
    
    # unpack data
    S11, S12, S13, S22, S23 = svec
    S33 = -S11 - S22
    
    # create matrix
    S = Matrix{Float32}(undef,3,3)
    
    # assign matrix entries
    S[1,1] = S11
    S[1,2] = S12
    S[2,1] = S12
    S[1,3] = S13
    S[3,1] = S13
    S[2,2] = S22
    S[2,3] = S23
    S[3,2] = S23
    S[3,3] = S33
    
    # normalize and return
    #return S ./ sqrt(sum(S.^2))
    return S ./ norm(S)    
    
end

### Given Stress Vector from Inversion, Assemble 3x3 Symmetric Tensor
#   Note, we are following Beauce's codes so s = [S11,S12,S13,S22,S23]
function assemble_stress_tensor(svec::Vector{Float64})
    
    # unpack data
    S11, S12, S13, S22, S23 = svec
    S33 = -S11 - S22
    
    # create matrix
    S = Matrix{Float64}(undef,3,3)
    
    # assign matrix entries
    S[1,1] = S11
    S[1,2] = S12
    S[2,1] = S12
    S[1,3] = S13
    S[3,1] = S13
    S[2,2] = S22
    S[2,3] = S23
    S[3,2] = S23
    S[3,3] = S33
    
    # normalize and return
    #return S ./ sqrt(sum(S.^2))
    return S ./ norm(S)
end

### Compute Principal Stresses
# Note that extension is positive, so we can keep the default sorting
function principal_stresses(S::Matrix{Float32})
    
    # eigenvector decomposition, returns in sorted order (negative to positive)
    E = eigen(S)
    
    # compute shape R (note extension is positive)
    R = (E.values[1]-E.values[2])/(E.values[1]-E.values[3])
    
    # return eigenvectors as columns, eigenvalues (extension is positive)
    return E.vectors, E.values, R
end

### Compute Principal Stresses
# Note that extension is positive, so we can keep the default sorting
function principal_stresses(S::Matrix{Float64})
    
    # eigenvector decomposition, returns in sorted order (negative to positive)
    E = eigen(S)
    
    # compute shape R (note extension is positive)
    R = (E.values[1]-E.values[2])/(E.values[1]-E.values[3])
    
    # return eigenvectors as columns, eigenvalues (extension is positive)
    return E.vectors, E.values, R
end

### Compute Azimuth and Plunge of a Vector (x1=N,x2=W,x3=U)
#  Azimuth is Angle between the north and the line.
#  Plunge is Angle between the horizontal plane and the line.
#  Calculations are for the lower hemisphere
function get_azi_plunge(v::Vector{Float32})
    
    # get correct sign the plunges into lower hemisphere
    if v[3] > 0.0f0
        u = -v
    else
        u = v
    end
    
    # in this coordinate system, note that the the trigonometric sense
    #  is opposite to the azimuthal sense, hence the -1
    azimuth = -1.0f0 * atand(u[2],u[1]) + 360.0f0
    
    # the plunge is measure downward from the end of the line of this bearing
    plunge = acosd(u[3]) - 90.0f0
    
    # return
    return azimuth % 360.0f0, plunge

end

### Compute Azimuth and Plunge of a Vector (x1=N,x2=W,x3=U)
#  Azimuth is Angle between the north and the line.
#  Plunge is Angle between the horizontal plane and the line.
#  Calculations are for the lower hemisphere
function get_azi_plunge(v::Vector{Float64})
    
    # get correct sign the plunges into lower hemisphere
    if v[3] > 0.0
        u = -v
    else
        u = v
    end
    
    # in this coordinate system, note that the the trigonometric sense
    #  is opposite to the azimuthal sense, hence the -1
    azimuth = -1.0 * atand(u[2],u[1]) + 360.0
    
    # the plunge is measure downward from the end of the line of this bearing
    plunge = acosd(u[3]) - 90.0
    
    # return
    return azimuth % 360.0, plunge

end

### Compute Shear and Normal Tractions, Given a Stress Tensor and Unit Normal Vector 
function normal_shear_tractions(S::Matrix{Float32},unorm::Matrix{Float32})
    
    # traction vector
    vecTR = unorm * S

    # normal traction vector
    normTR =  sum(vecTR .* unorm, dims=2) .* unorm
    
    # shear traction vector
    shearTR = vecTR .- normTR
   
    # return
    return normTR, shearTR
end

### Compute Shear and Normal Tractions, Given a Stress Tensor and Unit Normal Vector 
function normal_shear_tractions(S::Matrix{Float64},unorm::Matrix{Float64})
    
    # traction vector
    vecTR = unorm * S

    # normal traction vector
    normTR =  sum(vecTR .* unorm, dims=2) .* unorm
    
    # shear traction vector
    shearTR = vecTR .- normTR
   
    # return
    return normTR, shearTR
end

### Compute Misfit Between Shear Stress and Slip Direction
function angular_misfit(S::Matrix{Float32},strikes::Vector{Float32},
        dips::Vector{Float32},rakes::Vector{Float32})
    
    # compute the normal and slip vectors for plane
    xnorm = hcat(get_normal(strikes,dips)...) # unpacked
    xslip = hcat(get_slip(strikes,dips,rakes)...) # unpacked
    
    # compute normal and shear tractions
    _, shearTR = normal_shear_tractions(S,xnorm)
    
    # normalize the shear traction to a unit vector
    shearTR ./= sqrt.(sum(shearTR.^2,dims=2))
    
    # compute angle in degrees between slip and shear
    misfit = round_acosd(sum(xslip .* shearTR, dims=2))
    
    # return as a vector
    return vec(misfit) 
end
# --- similar to the above but returns just the mean
function mean_angular_misfit(S::Matrix{Float32},strikes::Vector{Float32},
    dips::Vector{Float32},rakes::Vector{Float32})

    # compute the normal and slip vectors for plane
    xnorm = hcat(get_normal(strikes,dips)...) # unpacked
    xslip = hcat(get_slip(strikes,dips,rakes)...) # unpacked

    # compute normal and shear tractions
    _, shearTR = normal_shear_tractions(S,xnorm)

    # normalize the shear traction to a unit vector
    shearTR ./= sqrt.(sum(shearTR.^2,dims=2))

    # compute angle in degrees between slip and shear
    return mean(round_acosd(sum(xslip .* shearTR, dims=2)))

end

### Compute Misfit Between Shear Stress and Slip Direction
function angular_misfit(S::Matrix{Float64},strikes::Vector{Float64},
    dips::Vector{Float64},rakes::Vector{Float64})

    # compute the normal and slip vectors for plane
    xnorm = hcat(get_normal(strikes,dips)...) # unpacked
    xslip = hcat(get_slip(strikes,dips,rakes)...) # unpacked

    # compute normal and shear tractions
    _, shearTR = normal_shear_tractions(S,xnorm)

    # normalize the shear traction to a unit vector
    shearTR ./= sqrt.(sum(shearTR.^2,dims=2))

    # compute angle in degrees between slip and shear
    misfit = round_acosd(sum(xslip .* shearTR, dims=2))

    # return as a vector
    return vec(misfit) 
end
# --- similar to the above but returns just the mean
function mean_angular_misfit(S::Matrix{Float64},strikes::Vector{Float64},
    dips::Vector{Float64},rakes::Vector{Float64})

    # compute the normal and slip vectors for plane
    xnorm = hcat(get_normal(strikes,dips)...) # unpacked
    xslip = hcat(get_slip(strikes,dips,rakes)...) # unpacked

    # compute normal and shear tractions
    _, shearTR = normal_shear_tractions(S,xnorm)

    # normalize the shear traction to a unit vector
    shearTR ./= sqrt.(sum(shearTR.^2,dims=2))

    # compute angle in degrees between slip and shear
    return mean(round_acosd(sum(xslip .* shearTR, dims=2)))

end

### Compute Instability
function compute_instability(V::Matrix{Float32},R::Float32, friction::Float32,
        strike1::Vector{Float32},dip1::Vector{Float32},rake1::Vector{Float32},
        strike2::Vector{Float32},dip2::Vector{Float32},rake2::Vector{Float32})
    
    # compute the normal and slip vectors for plane 1
    norm1 = hcat(get_normal(strike1,dip1)...) # unpacked
    slip1 = hcat(get_slip(strike1,dip1,rake1)...) # unpacked
 
    # compute the normal and slip vectors for plane 2
    norm2 = hcat(get_normal(strike2,dip2)...) # unpacked
    slip2 = hcat(get_slip(strike2,dip2,rake2)...) # unpacked
    
    # project into principal directions
    norm1 = norm1 * V
    norm2 = norm2 * V
    slip1 = slip1 * V
    slip2 = slip2 * V
    
    # normalized principal stresses - extension is positive!
    sig1 = -1.0f0
    sig2 = 2.0f0*R - 1.0f0
    sig3 = 1.0f0
    
    # compute normal stress on planes
    sigN1 = sig1.*(norm1[:,1].^2) + sig2.*(norm1[:,2].^2) + sig3.*(norm1[:,3].^2) 
    sigN2 = sig1.*(norm2[:,1].^2) + sig2.*(norm2[:,2].^2) + sig3.*(norm2[:,3].^2)
    
    # compute shear stress on planes
    sigS1 = round_sqrt((sig1^2).*(norm1[:,1].^2) .+ (sig2^2).*(norm1[:,2].^2) .+ (
            sig3^2).*(norm1[:,3].^2) .- sigN1.^2)
    sigS2 = round_sqrt((sig1^2).*(norm2[:,1].^2) .+ (sig2^2).*(norm2[:,2].^2) .+ (
            sig3^2).*(norm2[:,3].^2) .- sigN2.^2)
    
    # the critical values depend on friction
    tauC = 1.0f0 / sqrt(1.0f0 + friction^2)
    sigC = friction / sqrt(1.0f0 + friction^2)
    
    # compute instabilities
    Ic = tauC - friction * (sig1 - sigC) 
    I1 = (sigS1 .- friction.*(sig1 .- sigN1) ) / Ic
    I2 = (sigS2 .- friction.*(sig1 .- sigN2) ) / Ic
    
    # keep the most unstable planes
    ix1 = I1 .>= I2
    I = ifelse.(ix1,I1,I2)
    strike = ifelse.(ix1,strike1,strike2)
    dip = ifelse.(ix1,dip1,dip2)
    rake = ifelse.(ix1,rake1,rake2)
    
    # return
    return I, strike, dip, rake
    
end

### Compute Instability
function compute_instability(V::Matrix{Float64},R::Float64, friction::Float64,
    strike1::Vector{Float64},dip1::Vector{Float64},rake1::Vector{Float64},
    strike2::Vector{Float64},dip2::Vector{Float64},rake2::Vector{Float64})

    # compute the normal and slip vectors for plane 1
    norm1 = hcat(get_normal(strike1,dip1)...) # unpacked
    slip1 = hcat(get_slip(strike1,dip1,rake1)...) # unpacked

    # compute the normal and slip vectors for plane 2
    norm2 = hcat(get_normal(strike2,dip2)...) # unpacked
    slip2 = hcat(get_slip(strike2,dip2,rake2)...) # unpacked

    # project into principal directions
    norm1 = norm1 * V
    norm2 = norm2 * V
    slip1 = slip1 * V
    slip2 = slip2 * V

    # normalized principal stresses - extension is positive!
    sig1 = -1.0
    sig2 = 2.0*R - 1.0
    sig3 = 1.0

    # compute normal stress on planes
    sigN1 = sig1.*(norm1[:,1].^2) + sig2.*(norm1[:,2].^2) + sig3.*(norm1[:,3].^2) 
    sigN2 = sig1.*(norm2[:,1].^2) + sig2.*(norm2[:,2].^2) + sig3.*(norm2[:,3].^2)

    # compute shear stress on planes
    sigS1 = round_sqrt((sig1^2).*(norm1[:,1].^2) .+ (sig2^2).*(norm1[:,2].^2) .+ (
            sig3^2).*(norm1[:,3].^2) .- sigN1.^2)
    sigS2 = round_sqrt((sig1^2).*(norm2[:,1].^2) .+ (sig2^2).*(norm2[:,2].^2) .+ (
            sig3^2).*(norm2[:,3].^2) .- sigN2.^2)

    # the critical values depend on friction
    tauC = 1.0 / sqrt(1.0 + friction^2)
    sigC = friction / sqrt(1.0 + friction^2)

    # compute instabilities
    Ic = tauC - friction * (sig1 - sigC) 
    I1 = (sigS1 .- friction.*(sig1 .- sigN1) ) / Ic
    I2 = (sigS2 .- friction.*(sig1 .- sigN2) ) / Ic

    # keep the most unstable planes
    ix1 = I1 .>= I2
    I = ifelse.(ix1,I1,I2)
    strike = ifelse.(ix1,strike1,strike2)
    dip = ifelse.(ix1,dip1,dip2)
    rake = ifelse.(ix1,rake1,rake2)

    # return
    return I, strike, dip, rake

end

#### Similar to the Above but with Mechanism Matrix (N x 6: s1, d1, r1, s2, d2, r2)
function compute_instability(
    V::Matrix{Float32},R::Float32,friction::Float32,mechM::Matrix{Float32})

    # compute the normal and slip vectors for plane 1
    norm1 = hcat(get_normal(mechM[:,1],mechM[:,2])...) # unpacked
    slip1 = hcat(get_slip(mechM[:,1],mechM[:,2],mechM[:,3])...) # unpacked

    # compute the normal and slip vectors for plane 2
    norm2 = hcat(get_normal(mechM[:,4],mechM[:,5])...) # unpacked
    slip2 = hcat(get_slip(mechM[:,4],mechM[:,5],mechM[:,6])...) # unpacked

    # project into principal directions
    norm1 = norm1 * V
    norm2 = norm2 * V
    slip1 = slip1 * V
    slip2 = slip2 * V

    # normalized principal stresses - extension is positive!
    sig1 = -1.0f0
    sig2 = 2.0f0*R - 1.0f0
    sig3 = 1.0f0

    # compute normal stress on planes
    sigN1 = sig1.*(norm1[:,1].^2) + sig2.*(norm1[:,2].^2) + sig3.*(norm1[:,3].^2) 
    sigN2 = sig1.*(norm2[:,1].^2) + sig2.*(norm2[:,2].^2) + sig3.*(norm2[:,3].^2)

    # compute shear stress on planes
    sigS1 = round_sqrt((sig1^2).*(norm1[:,1].^2) .+ (sig2^2).*(norm1[:,2].^2) .+ (
            sig3^2).*(norm1[:,3].^2) .- sigN1.^2)
    sigS2 = round_sqrt((sig1^2).*(norm2[:,1].^2) .+ (sig2^2).*(norm2[:,2].^2) .+ (
            sig3^2).*(norm2[:,3].^2) .- sigN2.^2)

    # the critical values depend on friction
    tauC = 1.0f0 / sqrt(1.0f0 + friction^2)
    sigC = friction / sqrt(1.0f0 + friction^2)

    # compute instabilities
    Ic = tauC - friction * (sig1 - sigC) 
    I1 = (sigS1 .- friction.*(sig1 .- sigN1) ) / Ic
    I2 = (sigS2 .- friction.*(sig1 .- sigN2) ) / Ic

    # keep the most unstable planes
    ix1 = I1 .>= I2
    I = ifelse.(ix1,I1,I2)
    strike = ifelse.(ix1,mechM[:,1],mechM[:,4])
    dip = ifelse.(ix1,mechM[:,2],mechM[:,5])
    rake = ifelse.(ix1,mechM[:,3],mechM[:,6])

    # return
    return I, strike, dip, rake

end

function compute_instability(
    V::Matrix{Float64},R::Float64, friction::Float64,mechM::Matrix{Float64})

    # compute the normal and slip vectors for plane 1
    norm1 = hcat(get_normal(mechM[:,1],mechM[:,2])...) # unpacked
    slip1 = hcat(get_slip(mechM[:,1],mechM[:,2],mechM[:,3])...) # unpacked

    # compute the normal and slip vectors for plane 2
    norm2 = hcat(get_normal(mechM[:,4],mechM[:,5])...) # unpacked
    slip2 = hcat(get_slip(mechM[:,4],mechM[:,5],mechM[:,6])...) # unpacked

    # project into principal directions
    norm1 = norm1 * V
    norm2 = norm2 * V
    slip1 = slip1 * V
    slip2 = slip2 * V

    # normalized principal stresses - extension is positive!
    sig1 = -1.0
    sig2 = 2.0*R - 1.0
    sig3 = 1.0

    # compute normal stress on planes
    sigN1 = sig1.*(norm1[:,1].^2) + sig2.*(norm1[:,2].^2) + sig3.*(norm1[:,3].^2) 
    sigN2 = sig1.*(norm2[:,1].^2) + sig2.*(norm2[:,2].^2) + sig3.*(norm2[:,3].^2)

    # compute shear stress on planes
    sigS1 = round_sqrt((sig1^2).*(norm1[:,1].^2) .+ (sig2^2).*(norm1[:,2].^2) .+ (
            sig3^2).*(norm1[:,3].^2) .- sigN1.^2)
    sigS2 = round_sqrt((sig1^2).*(norm2[:,1].^2) .+ (sig2^2).*(norm2[:,2].^2) .+ (
            sig3^2).*(norm2[:,3].^2) .- sigN2.^2)

    # the critical values depend on friction
    tauC = 1.0 / sqrt(1.0 + friction^2)
    sigC = friction / sqrt(1.0 + friction^2)

    # compute instabilities
    Ic = tauC - friction * (sig1 - sigC) 
    I1 = (sigS1 .- friction.*(sig1 .- sigN1) ) / Ic
    I2 = (sigS2 .- friction.*(sig1 .- sigN2) ) / Ic

    # keep the most unstable planes
    ix1 = I1 .>= I2
    I = ifelse.(ix1,I1,I2)
    strike = ifelse.(ix1,mechM[:,1],mechM[:,4])
    dip = ifelse.(ix1,mechM[:,2],mechM[:,5])
    rake = ifelse.(ix1,mechM[:,3],mechM[:,6])

    # return
    return I, strike, dip, rake

end

### Rounding functions for acos and sqrt ###

# keeps argument positive
function round_sqrt(x::Vector{Float32})
    return sqrt.(ifelse.(x.>=0.0f0,x,0.0f0))
end

# keeps argument positive
function round_sqrt(x::Vector{Float64})
    return sqrt.(ifelse.(x.>=0.0,x,0.0))
end

# keeps argument positive
function round_sqrt(x::Matrix{Float32})
    return sqrt.(ifelse.(x.>=0.0f0,x,0.0f0))
end

# keeps argument positive
function round_sqrt(x::Matrix{Float64})
    return sqrt.(ifelse.(x.>=0.0,x,0.0))
end

# keeps argument between [-1,1]
function round_acosd(x::Vector{Float32})
    return acosd.(ifelse.(abs.(x).>=1.0f0,sign.(x),x))
end

# keeps argument between [-1,1]
function round_acosd(x::Vector{Float64})
    return acosd.(ifelse.(abs.(x).>=1.0,sign.(x),x))
end

# keeps argument between [-1,1]
function round_acosd(x::Matrix{Float32})
    return acosd.(ifelse.(abs.(x).>=1.0f0,sign.(x),x))
end

# keeps argument between [-1,1]
function round_acosd(x::Matrix{Float64})
    return acosd.(ifelse.(abs.(x).>=1.0,sign.(x),x))
end


#### mt2sdr - float 32
#     Adapted from obspy.imaging.beachball, originally from MATLAB script
# `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
# written by Andy Michael, Chen Ji and Oliver Boyd.
function mt2sdr(mrr::Float32,mtt::Float32,mpp::Float32,mrt::Float32,mrp::Float32,mtp::Float32)

    # make MT matrix
    MT = [[mrr mrt mrp]
          [mrt mtt mtp]
          [mrp mtp mpp]]

    # eigenvector decomposition, returns in sorted order (negative to positive)
    E = eigen(MT) # E.values is a vector, E.vectors is matrix 

    # new coordinate system
    d = [E.values[2],E.values[1],E.values[3]]
    v = [[ E.vectors[2, 2] -E.vectors[2, 1] -E.vectors[2, 3]]
         [ E.vectors[3, 2] -E.vectors[3, 1] -E.vectors[3, 3]]
         [-E.vectors[1, 2]  E.vectors[1, 1]  E.vectors[1, 3]]]

    # find min / max positions
    imin, imax = argmin(d), argmax(d)

    # define vectors
    ae = (v[1:3, imax] + v[1:3, imin]) / sqrt(2.0f0)
    an = (v[1:3, imax] - v[1:3, imin]) / sqrt(2.0f0)

    # normalization
    aer = sqrt(ae[1]^2 + ae[2]^2 + ae[3]^2)
    anr = sqrt(an[1]^2 + an[2]^2 + an[3]^2)
    ae = ae / aer
    an = an / anr

    # sign convention
    if an[3] <= 0.0f0
        an1 = an
        ae1 = ae
    else
        an1 = -an
        ae1 = -ae
    end

    # calculate nodal planes
    ft, fd, fl = tdl(an1,ae1)
    strike1, dip1, rake1 = 360.0f0 - ft, fd, 180.0f0 - fl
    strike2, dip2, rake2 = aux_plane([strike1],[dip1],[rake1])

    # return
    return strike1, dip1, rake1, strike2[1], dip2[1], rake2[1]

end

#### mt2sdr - float 64
#     Adapted from obspy.imaging.beachball, originally from MATLAB script
# `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
# written by Andy Michael, Chen Ji and Oliver Boyd.
function mt2sdr(mrr::Float64,mtt::Float64,mpp::Float64,mrt::Float64,mrp::Float64,mtp::Float64)

    # make MT matrix
    MT = [[mrr mrt mrp]
          [mrt mtt mtp]
          [mrp mtp mpp]]

    # eigenvector decomposition, returns in sorted order (negative to positive)
    E = eigen(MT) # E.values is a vector, E.vectors is matrix 

    # new coordinate system
    d = [E.values[2],E.values[1],E.values[3]]
    v = [[ E.vectors[2, 2] -E.vectors[2, 1] -E.vectors[2, 3]]
         [ E.vectors[3, 2] -E.vectors[3, 1] -E.vectors[3, 3]]
         [-E.vectors[1, 2]  E.vectors[1, 1]  E.vectors[1, 3]]]

    # find min / max positions
    imin, imax = argmin(d), argmax(d)

    # define vectors
    ae = (v[1:3, imax] + v[1:3, imin]) / sqrt(2.0)
    an = (v[1:3, imax] - v[1:3, imin]) / sqrt(2.0)

    # normalization
    aer = sqrt(ae[1]^2 + ae[2]^2 + ae[3]^2)
    anr = sqrt(an[1]^2 + an[2]^2 + an[3]^2)
    ae = ae / aer
    an = an / anr

    # sign convention
    if an[3] <= 0.0
        an1 = an
        ae1 = ae
    else
        an1 = -an
        ae1 = -ae
    end

    # calculate nodal planes
    ft, fd, fl = tdl(an1,ae1)
    strike1, dip1, rake1 = 360.0 - ft, fd, 180.0 - fl
    strike2, dip2, rake2 = aux_plane([strike1],[dip1],[rake1])

    # return
    return strike1, dip1, rake1, strike2[1], dip2[1], rake2[1]

end

### vectorized versions
function mt2sdr(mrr::Vector{Float32},mtt::Vector{Float32},mpp::Vector{Float32},
    mrt::Vector{Float32},mrp::Vector{Float32},mtp::Vector{Float32})

    # allocate arrays
    N = length(mrr)
    strike1, dip1, rake1 = Vector{Float32}(undef,N), Vector{Float32}(undef,N), Vector{Float32}(undef,N)
    strike2, dip2, rake2 = Vector{Float32}(undef,N), Vector{Float32}(undef,N), Vector{Float32}(undef,N)

    # calculations for each entry
    for ii in 1:N
        strike1[ii], dip1[ii], rake1[ii], strike2[ii], dip2[ii], rake2[ii] = mt2sdr(
            mrr[ii], mtt[ii], mpp[ii], mrt[ii], mrp[ii], mtp[ii])
    end

    # return
    return strike1, dip1, rake1, strike2, dip2, rake2
end
function mt2sdr(mrr::Vector{Float64},mtt::Vector{Float64},mpp::Vector{Float64},
    mrt::Vector{Float64},mrp::Vector{Float64},mtp::Vector{Float64})

    # allocate arrays
    N = length(mrr)
    strike1, dip1, rake1 = Vector{Float64}(undef,N), Vector{Float64}(undef,N), Vector{Float64}(undef,N)
    strike2, dip2, rake2 = Vector{Float64}(undef,N), Vector{Float64}(undef,N), Vector{Float64}(undef,N)

    # calculations for each entry
    for ii in 1:N
        strike1[ii], dip1[ii], rake1[ii], strike2[ii], dip2[ii], rake2[ii] = mt2sdr(
            mrr[ii], mtt[ii], mpp[ii], mrt[ii], mrp[ii], mtp[ii])
    end

    # return
    return strike1, dip1, rake1, strike2, dip2, rake2
end

#### helper function for mt2sdr - Float32
#     Adapted from obspy.imaging.beachball, originally from MATLAB script
# `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
# written by Andy Michael, Chen Ji and Oliver Boyd.
function tdl(an::Vector{Float32}, bn::Vector{Float32})
    """
    Helper function for mt2plane.

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.
    """

    # setup
    xn = an[1]
    yn = an[2]
    zn = an[3]
    xe = bn[1]
    ye = bn[2]
    ze = bn[3]
    aaa = 1.0f0 / (1000000.0f0)
    con = 57.2957795f0

    # calculations
    if abs(zn) < aaa
        fd = 90.0f0
        axn = abs(xn)
        if axn>1.0f0; axn = 1.0f0; end
        ft = asin(axn) * con
        st = -xn
        ct = yn
        if (st >= 0.0f0) & (ct < 0.0f0); ft = 180.0f0 - ft; end
        if (st < 0.0f0) & (ct <= 0.0f0); ft = 180.0f0 + ft; end
        if (st < 0.0f0) & (ct > 0.0f0); ft = 360.0f0 - ft; end
        fl = asin(abs(ze)) * con
        sl = -ze
        if abs(xn) < aaa
            cl = xe / yn
        else
            cl = -ye / xn
        end
        if (sl >= 0.0f0) & (cl < 0.0f0); fl = 180.0f0 - fl; end
        if (sl < 0.0f0) & (cl <= 0.0f0); fl = fl - 180.0f0; end
        if (sl < 0.0f0) & (cl > 0.0f0); fl = -fl; end
    else
        if -zn > 1.0f0; zn = -1.0f0; end
        fdh = acos(-zn)
        fd = fdh * con
        sd = sin(fdh)
        if sd == 0.0f0; return; end
        st = -xn / sd
        ct = yn / sd
        sx = abs(st)
        if sx > 1.0f0; sx = 1.0f0; end
        ft = asin(sx) * con
        if (st >= 0.0f0) & (ct < 0.0f0); ft = 180.0f0 - ft; end
        if (st < 0.0f0) & (ct <= 0.0f0); ft = 180.0f0 + ft; end
        if (st < 0.0f0) & (ct > 0.0f0); ft = 360.0f0 - ft; end
        sl = -ze / sd
        sx = abs(sl)
        if sx > 1.0f0; sx = 1.0f0; end
        fl = asin(sx) * con
        if st == 0.0f0
            cl = xe / ct
        else
            xxx = yn * zn * ze / sd / sd + ye
            cl = -sd * xxx / xn
            if ct == 0.0f0; cl = ye / st; end
        end
        if (sl >= 0.0f0) & (cl < 0.0f0); fl = 180.0f0 - fl; end
        if (sl < 0.0f0) & (cl <= 0.0f0);  fl = fl - 180.0f0; end
        if (sl < 0.0f0) & (cl > 0.0f0); fl = -fl; end
    end
    
    # return
    return ft, fd, fl

end

#### helper function for mt2sdr - Float64
#     Adapted from obspy.imaging.beachball, originally from MATLAB script
# `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
# written by Andy Michael, Chen Ji and Oliver Boyd.
function tdl(an::Vector{Float64}, bn::Vector{Float64})
    """
    Helper function for mt2plane.

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.
    """

    # setup
    xn = an[1]
    yn = an[2]
    zn = an[3]
    xe = bn[1]
    ye = bn[2]
    ze = bn[3]
    aaa = 1.0 / (1000000.0)
    con = 57.2957795

    # calculations
    if abs(zn) < aaa
        fd = 90.0
        axn = abs(xn)
        if axn>1.0; axn = 1.0; end
        ft = asin(axn) * con
        st = -xn
        ct = yn
        if (st >= 0.0) & (ct < 0.0); ft = 180.0 - ft; end
        if (st < 0.0) & (ct <= 0.0); ft = 180.0 + ft; end
        if (st < 0.0) & (ct > 0.0); ft = 360.0 - ft; end
        fl = asin(abs(ze)) * con
        sl = -ze
        if abs(xn) < aaa
            cl = xe / yn
        else
            cl = -ye / xn
        end
        if (sl >= 0.0) & (cl < 0.0); fl = 180.0 - fl; end
        if (sl < 0.0) & (cl <= 0.0); fl = fl - 180.0; end
        if (sl < 0.0) & (cl > 0.0); fl = -fl; end
    else
        if -zn > 1.0; zn = -1.0; end
        fdh = acos(-zn)
        fd = fdh * con
        sd = sin(fdh)
        if sd == 0.0; return; end
        st = -xn / sd
        ct = yn / sd
        sx = abs(st)
        if sx > 1.0; sx = 1.0; end
        ft = asin(sx) * con
        if (st >= 0.0) & (ct < 0.0); ft = 180.0 - ft; end
        if (st < 0.0) & (ct <= 0.0); ft = 180.0 + ft; end
        if (st < 0.0) & (ct > 0.0); ft = 360.0 - ft; end
        sl = -ze / sd
        sx = abs(sl)
        if sx > 1.0; sx = 1.0; end
        fl = asin(sx) * con
        if st == 0.0
            cl = xe / ct
        else
            xxx = yn * zn * ze / sd / sd + ye
            cl = -sd * xxx / xn
            if ct == 0.0; cl = ye / st; end
        end
        if (sl >= 0.0) & (cl < 0.0); fl = 180.0 - fl; end
        if (sl < 0.0) & (cl <= 0.0);  fl = fl - 180.0; end
        if (sl < 0.0) & (cl > 0.0); fl = -fl; end
    end
    
    # return
    return ft, fd, fl

end

### Functions to convert strike, dip, rake --> MT in spherical coordinates
# - Adapted from MomentTensors.jl
function sdr2mt(s::Float32, d::Float32, r::Float32, M0::Float32)
    tt = -M0*(sind(d)*cosd(r)*sind(2s) + sind(2d)*sind(r)*sind(s)^2)
    tp = -M0*(sind(d)*cosd(r)*cosd(2s) + sind(2d)*sind(r)*sind(2s)/2)
    rt = -M0*(cosd(d)*cosd(r)*cosd(s)  + cosd(2d)*sind(r)*sind(s))
    pp =  M0*(sind(d)*cosd(r)*sind(2s) - sind(2d)*sind(r)*cosd(s)^2)
    rp =  M0*(cosd(d)*cosd(r)*sind(s)  - cosd(2d)*sind(r)*cosd(s))
    rr =  M0*sind(2d)*sind(r)
    return rr, tt, pp, rt, rp, tp
end
function sdr2mt(s::Float64, d::Float64, r::Float64, M0::Float64)
    tt = -M0*(sind(d)*cosd(r)*sind(2s) + sind(2d)*sind(r)*sind(s)^2)
    tp = -M0*(sind(d)*cosd(r)*cosd(2s) + sind(2d)*sind(r)*sind(2s)/2)
    rt = -M0*(cosd(d)*cosd(r)*cosd(s)  + cosd(2d)*sind(r)*sind(s))
    pp =  M0*(sind(d)*cosd(r)*sind(2s) - sind(2d)*sind(r)*cosd(s)^2)
    rp =  M0*(cosd(d)*cosd(r)*sind(s)  - cosd(2d)*sind(r)*cosd(s))
    rr =  M0*sind(2d)*sind(r)
    return rr, tt, pp, rt, rp, tp
end
function sdr2mt(s::Vector{Float32}, d::Vector{Float32}, r::Vector{Float32}, M0::Vector{Float32})
    tt = -M0.*(sind.(d).*cosd.(r).*sind.(2s) .+ sind.(2d).*sind.(r)*sind.(s).^2)
    tp = -M0.*(sind.(d).*cosd.(r).*cosd.(2s) .+ sind.(2d).*sind.(r)*sind.(2s)./2)
    rt = -M0.*(cosd.(d).*cosd.(r).*cosd.(s)  .+ cosd.(2d).*sind.(r)*sind.(s))
    pp =  M0.*(sind.(d).*cosd.(r).*sind.(2s) .- sind.(2d).*sind.(r)*cosd.(s).^2)
    rp =  M0.*(cosd.(d).*cosd.(r).*sind.(s)  .- cosd.(2d).*sind.(r)*cosd.(s))
    rr =  M0.*sind.(2d).*sind.(r)
    return rr, tt, pp, rt, rp, tp
end
function sdr2mt(s::Vector{Float64}, d::Vector{Float64}, r::Vector{Float64}, M0::Vector{Float64})
    tt = -M0.*(sind.(d).*cosd.(r).*sind.(2s) .+ sind.(2d).*sind.(r)*sind.(s).^2)
    tp = -M0.*(sind.(d).*cosd.(r).*cosd.(2s) .+ sind.(2d).*sind.(r)*sind.(2s)./2)
    rt = -M0.*(cosd.(d).*cosd.(r).*cosd.(s)  .+ cosd.(2d).*sind.(r)*sind.(s))
    pp =  M0.*(sind.(d).*cosd.(r).*sind.(2s) .- sind.(2d).*sind.(r)*cosd.(s).^2)
    rp =  M0.*(cosd.(d).*cosd.(r).*sind.(s)  .- cosd.(2d).*sind.(r)*cosd.(s))
    rr =  M0.*sind.(2d).*sind.(r)
    return rr, tt, pp, rt, rp, tp
end
