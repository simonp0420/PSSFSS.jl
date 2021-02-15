module GSMs
export GSM, cascade, cascade!, initialize_gsm_file, append_gsm_data, read_gsm_file

using LinearAlgebra
using StaticArrays: SA
using OffsetArrays
using Unitful: @u_str, ustrip, unit
using JLD2
using ..Substrate: Layer, TEorTM, TE, TM, Gblock
using ..PGF: mysqrt
using ..Constants: η₀, min_elength
using ..Sheets: Sheet, RWGSheet, find_unique_periods
using FileIO: load
  
mutable struct GSM{T1 <: AbstractMatrix, T2 <: AbstractMatrix}
    s11::T1
    s12::T2
    s21::T2
    s22::T1
end

function GSM(n1::Int, n2::Int)
    gsm = GSM(zeros(ComplexF64, n1, n1), zeros(ComplexF64, n1, n2),
              zeros(ComplexF64, n2, n1), zeros(ComplexF64, n2, n2))
    gsm.s12[diagind(gsm.s12)] .= one(eltype(gsm.s12))
    gsm.s21[diagind(gsm.s21)] .= one(eltype(gsm.s21))
    gsm
end


"""
    cascade(a::GSM, b::GSM) -> c::GSM

Cascade a pair of generalized scattering matrices (GSMs).

## Input arguments

- `a`, `b`:  variables containing the input GSMs that are to be cascaded.  The matrices 
must be conformable, i.e., `n2a == n1b`, where `n2a` is the number of modes in Region
2 for GSM `a`, and `n1b` is the number of modes in Region 1 of GSM `b`.
"""
function cascade(a::GSM, b::GSM)::GSM
    n2a = size(a.s22, 2)
    n1b = size(b.s11, 1)
    n1b ≠ n2a && error("Non-conformable arrays")
    n1 =  size(a.s11, 1)
    n2 =  size(b.s22, 2)
    ninner = n1b

    # Equation (3.35) of the theory documentation:
    gprod1 = (I - a.s22 * b.s11) \ a.s21
    s21 = b.s21 * gprod1
    s11 = a.s11 + (a.s12 * (b.s11 * gprod1))
    gprod2 = (I - b.s11 * a.s22) \ b.s12
    s12 = a.s12 * gprod2
    s22 = b.s22 + (b.s21 * (a.s22 * gprod2))
    GSM(s11, s12, s21, s22)
end



"""
    cascade!(a::GSM, layer::Layer, t=NaN)

Cascade a GSM (generalized scattering matrix) on the left with a 
dielectric slab on the right. On exit, `a` is modified. `t` is an
layer thickness in layer units.  If not passed or set to `NaN`, 
the thickness given in the `layer` variable is used.  It is assumed that
the modes corresponding to the side `2` ports of the GSM are already properly 
defined (i.e. normalized) consistent with the dielectric slab electrical properties.
"""
function cascade!(a::GSM, layer::Layer, t=NaN)
    n1 = size(a.s11, 1)
    n2 = size(a.s22, 1)
    n2 ≠ length(layer.γ) && @error "# modes not consistent" n2 length(layer.γ) exception=ErrorException
    # Obtain layer thickness in meters:
    if !isnan(t)
        units_per_meter = ustrip(Float64, unit(layer.user_width), 1u"m")
        d = t / units_per_meter        # Layer thickness in meters.
    else
        d = layer.width
    end
    d == 0 && return nothing
    #  Loop over each of the modes:
    for i in 1:n2
        p = exp(-layer.γ[i] * d)
        a.s12[:,i] .*= p
        a.s21[i,:] .*= p
        a.s22[:,i] .*= p
        a.s22[i,:] .*= p
    end
    return nothing
end


"z-component of the cross product of two 3-vectors"
zdotcross(a,b) = a[1] * b[2] - a[2] * b[1]



"""
    gsm_electric_gblock(layers::Vector{Layer}, s::Integer, k0::Float64) -> (gsm, tlgfvi, vincs)

Calculate the GSM (generalized scattering matrix) of a `GBLOCK` containing a single electric-type
FSS surface.  Also compute the quantities from Section 6.2 of the theory documentation needed to 
compute incident and scattered fields.

## Arguments

- `layers`:  Contains the layer parameters for the cascade structure.  Note that the first
 and last layer's thicknesses are not accounted for in this function.  They are assumed to be semi-infinite.
- `s`: Interface number (within layers) at which the FSS is located.
- `k0`: Free-space wavenumber in rad/m.

## Return Values

- `gsm::GSM`: A variable containing the generalized scattering matrix initialized by this function.

- `tlgfvi::Matrix{ComplexF64}` `tlgfvi` is the voltage transmission line Green's function 
 (TLGF) due to a unit current. `tlgfvi[q,1]` is the TLGF for mode `q` of Region 1, with source 
at ``z=z_s`` and observation point at ``z=z_1``. `tlgfvi[q,2]` is the TLGF for mode `q` of 
Region N, with source at ``z=z_s`` and observation point at ``z=z_{N-1}``. 

- `vincs::Matrix{ComplexF64}`: `vincs` is the transmission line voltage evaluated at ``z = z_s``
due to a Thevenin voltage source ``V_g = 2`` at either the left (``z=z_1``) or right (``z=z_{N-1}``)
ends of the equivalent circuit.   `vincs[q,1]` is the voltage for mode `q` with source at ``z=z_1``, 
`vincs[q,2]` is the voltage for mode `q` with the source at ``z=z_{N-1}``.

"""
function gsm_electric_gblock(layers::Vector{Layer}, s::Integer, k0::Float64)
    N = length(layers)
    cosha = OffsetArray(zeros(ComplexF64, N-2), 2:N-1)
    sinha = OffsetArray(zeros(ComplexF64, N-2), 2:N-1)
    tanha = OffsetArray(zeros(ComplexF64, N-2), 2:N-1)
    Z0 = zeros(ComplexF64, N)
    Zleft = zeros(ComplexF64, N)
    Zright = zeros(ComplexF64, N)
    
    k0sq = k0 * k0
    n1 = length(layers[1].P)  # Number of modes in Region 1.
    n2 = length(layers[N].P)  # Number of modes in Region N.
    gsm = GSM(n1,n2)
    n12max = max(n1, n2)
    vincs = zeros(ComplexF64, n12max, 2)
    tlgfvi = zeros(ComplexF64, n12max, 2)
    
    
    # Compute square root of the ratio of the Region 1 and Region n
    # unit cell areas.  This is needed to account for the fact that the 
    # two regions may not employ the same lattice to define their modes, 
    # hence the mode normalizations may not be compatible.
    root_area1oN = sqrt(abs(
        zdotcross(layers[N].β₁, layers[N].β₂) / zdotcross(layers[1].β₁, layers[1].β₂)))
    
    #  Loop over the region 1 source modes:
    for q1 in 1:n1
        Z0[1] = 1 / layers[1].Y[q1]  # Obtain region 1 modal impedance.
        Zleft[1] = Z0[1]  # Eq. (6.12a) of the theory documentation
        # Compute mag squared of transverse wavenumber (rad^2/m^2):
        β² = layers[1].β[q1] ⋅ layers[1].β[q1] 
        # Loop over each region:
        for i in 2:N
            # Obtain wavenumber (Region i) squared (rad^2/m^2):
            k² = k0sq * layers[i].μᵣ * layers[i].ϵᵣ
            γ = mysqrt(β² - k²) # Atten. constant (neper/m).
            
            if layers[1].P[q1] == TE
                Z0[i] = (im * (k0*η₀) * layers[i].μᵣ) / γ # TE
            else
                Z0[i] = γ / (im * (k0/η₀) * layers[i].ϵᵣ) # TM
            end
            if i < N
                # Compute remaining quantities for internal layers:
                cosha[i], sinha[i] = sincos(im * layers[i].width * γ)
                sinha[i] = -sinha[i]
                tanha[i] = sinha[i] / cosha[i]
                Zleft[i] = Z0[i] *
                    (Zleft[i-1] + Z0[i]*tanha[i]) / (Z0[i] + Zleft[i-1] * tanha[i])  # Eq. (6.12b).
            end
        end
        # Loop over regions in backwards order to get Zright:
        Zright[N-1] = Z0[N] # Eq. (6.4a)
        for i in N-2:-1:1
            Zright[i] = Z0[i+1] *
                (Zright[i+1] + Z0[i+1]*tanha[i+1]) / 
                (Z0[i+1] + Zright[i+1] * tanha[i+1])  # Eq. (6.4b).
        end
        
        # Compute voltage and current at z = z_1 using Eq. (6.3):
        current = 2  / (Zright[1] + Z0[1])
        voltage = current * Zright[1]
        # Compute reflection term using Eq. (6.2):
        gsm.s11[q1,q1] = voltage - 1
        
        # Step voltage and current up to junction s using Eq. (6.5):
        for i in 2:s
            voltage = voltage * cosha[i] - Z0[i] * current * sinha[i]
            current = voltage / Zright[i]
        end
        vincs[q1,1] = voltage  # Save incident voltage
        
        # Locate identical mode in region n:
        q2 = find_mode_index(layers[1].P[q1], layers[1].M[q1], layers[1].N[q1], layers[N])
        if q2 ≠ 0  # Do we need to compute S21 entry?
            # Step voltage and current up to z = z_{N-1} junction using Eq. (6.5):
            for i in s+1:N-1
                voltage = voltage * cosha[i] - Z0[i] * current * sinha[i]
                current = voltage / Zright[i]
            end
            # Compute transmission S-parameter using Eq. (8):
            gsm.s21[q2,q1] = voltage * layers[1].c[q1] / layers[N].c[q2] * root_area1oN
        end
        
        # Now compute TLGF's:
        current =  -Zright[s] / (Zleft[s] + Zright[s]) # Eqs. (6.19c) and (6.19a)
        voltage = -Zleft[s] * current                  # (6.19a)
        # Step down to z_1 using (6.13):
        for i in s-1:-1:1
            voltage = voltage * cosha[i+1] + Z0[i+1] * current * sinha[i+1]
            current = -voltage / Zleft[i]
        end
        tlgfvi[q1,1] = voltage # Save quantity used in (18a).
        
    end # Done with Region 1 modes.
    
    for q2 in 1:n2  #  Loop over the region N source modes:
        Z0[N] = 1 / layers[N].Y[q2]  # Obtain region N modal impedance.
        Zright[N-1] = Z0[N]  # Eq. (6.4a)
        # Compute mag squared of transverse wavenumber (rad^2/m^2):
        β² = layers[N].β[q2] ⋅ layers[N].β[q2]
        for i in N-1:-1:1 # Loop over each region
            k² = k0sq * layers[i].μᵣ * layers[i].ϵᵣ
            γ = mysqrt(β² - k²) # Atten. constant (neper/m).
            if layers[N].P[q2] == TE
                Z0[i] = (im * (k0*η₀) * layers[i].μᵣ) / γ
            else
                Z0[i] = γ / (im * (k0/η₀) * layers[i].ϵᵣ) 
            end
            if i > 1
                # Compute remaining quantities for internal layers:
                cosha[i], sinha[i] = sincos(im * layers[i].width * γ)
                sinha[i] = -sinha[i]
                tanha[i] = sinha[i] / cosha[i]
                Zright[i-1] = Z0[i] *
                    (Zright[i] + Z0[i]*tanha[i]) / (Z0[i] + Zright[i] * tanha[i])  # Eq. (6.4b)
            end
        end
        # Loop over regions in forward order to get Zleft:
        Zleft[1] = Z0[1] # Eq. (6.12a)
        for i in 2:N-1
            Zleft[i] = Z0[i] *
                (Zleft(i-1) + Z0[i]*tanha[i]) / (Z0[i] + Zleft(i-1) * tanha[i]) # Eq. (6.12b)
        end
        # Compute voltage and current at z = z_1 using Eq. (6.11):
        current = -2  / (Zleft[N-1] + Z0[N])
        voltage = -current * Zleft[N-1]
        gsm.s22[q2,q2] = voltage - 1 # Compute reflection term using Eq. (6.10)
        
        # Step voltage and current down to junction s using
        # Equation (6.13):
        for i in N-2:-1:s
            voltage = voltage * cosha[i+1] + Z0[i+1] * current * sinha[i+1]
            current = -voltage / Zleft[i]
        end
        vincs[q2,2] = voltage  # Save incident voltage at FSS plane.
        
        # Locate identical mode in region 1:
        q1 = find_mode_index(layers[N].P[q2], layers[N].M[q2], layers[N].N[q2], layers[1])
        if q1 ≠ 0
            # Step voltage and current down to z = z_{1} junction using Eq. (6.13):
            for i in s-1:-1:1
                voltage = voltage * cosha[i+1] + Z0[i+1] * current * sinha[i+1]
                current = -voltage / Zleft[i]
            end
            # Compute transmission S-parameter using Eq. (6.8)
            gsm.s12[q1,q2] = voltage * layers[N].c[q2] / layers[1].c[q1] / root_area1oN
        end
        
        # Now compute TLGF's:
        current =  Zleft[s] / (Zleft[s] + Zright[s]) # Eq. (6.19b)
        voltage = Zright[s] * current                  # Eq. (6.19a)
        # Step up to z_{N-1} using (6.13):
        for i in s+1:N-1
            voltage = voltage * cosha[i] - Z0[i] * current * sinha[i]
            current = voltage / Zright[i]
        end
        tlgfvi[q2,2] = voltage # Save quantity used in (6.18b)
    end  # Done with Region N modes
    return (gsm, tlgfvi, vincs)
end


"""
    gsm_magnetic_gblock(layers::Vector{Layer}, s::Integer, k0::Float64) -> (gsm, tlgfiv, iincs)

Calculate the GSM (generalized scattering matrix) of a `GBLOCK` containing a single magnetic-type
FSS surface.  Also compute the quantities from Section 6.3 of the theory documentation needed to 
compute incident and scattered fields.

## Arguments

- `layers`:  Contains the layer parameters for the cascade structure.  Note that the first
 and last layer's thicknesses are not accounted for in this function.  They are assumed to be semi-infinite.
- `s`: Interface number (within layers) at which the FSS is located.
- `k0`: Free-space wavenumber in rad/m.

## Return Values

- `gsm::GSM`: A variable containing the generalized scattering matrix initialized by this function.

- `tlgfiv::Matrix{ComplexF64}` `tlgfiv` is the current transmission line Green's function 
 (TLGF) due to a unit voltage. `tlgfiv[q,1]` is the TLGF for mode `q` of Region 1, with source 
at ``z=z_s`` and observation point at ``z=z_1``. `tlgfvi[q,2]` is the TLGF for mode `q` of 
Region N, with source at ``z=z_s`` and observation point at ``z=z_{N-1}``. 

- `iincs::Matrix{ComplexF64}`: `iincs` is the transmission line current evaluated at ``z = z_s``
due to a Norton current source ``I_g = 2`` at either the left (``z=z_1``) or right (``z=z_{N-1}``)
ends of the equivalent circuit.   `iincs[q,1]` is the current for mode `q` with source at ``z=z_1``, 
`iincs[q,2]` is the voltage for mode `q` with the source at ``z=z_{N-1}``.
"""
function gsm_magnetic_gblock(layers::Vector{Layer}, s::Integer, k0::Float64)
    N = length(layers)
    cosha = OffsetArray(zeros(ComplexF64, N-2), 2:N-1)
    sinha = OffsetArray(zeros(ComplexF64, N-2), 2:N-1)
    tanha = OffsetArray(zeros(ComplexF64, N-2), 2:N-1)
    Z0 = zeros(ComplexF64, N)
    Zleft = zeros(ComplexF64, N)
    Zright = zeros(ComplexF64, N)

    k0sq = k0 * k0
    n1 = length(layers[1].P)  # Number of modes in Region 1.
    n2 = length(layers[N].P)  # Number of modes in Region N.
    gsm = GSM(n1, n2)
    n12max = max(n1,n2)
    iincs = zeros(ComplexF64, n12max, 2)
    tlgfiv = zeros(ComplexF64, n12max, 2)

    # Loop over the region 1 source modes:
    for q1 in 1:n1
        Z0[1] = 1 / layers[1].Y[q1]  # Obtain region 1 modal impedance
        Zleft[1] = Z0[1]  # Semi-infinite line.
        # Compute mag squared of transverse wavenumber (rad^2/m^2):
        β² = layers[1].β[q1] ⋅ layers[1].β[q1] 
        # Loop over each region:
        for i in 2:s
            # Obtain wavenumber (Region i) squared (rad^2/m^2):
            k² = k0sq * layers[i].μᵣ * layers[i].ϵᵣ
            γ = mysqrt(β² - k²) # Atten. constant (neper/m).
            if layers[1].P[q1] == TE
                Z0[i] = (im * (k0*η₀) * layers[i].μᵣ) / γ # TE
            else
                Z0[i] = γ / (im * (k0/η₀) * layers[i].ϵᵣ) # TM
            end
            # Compute remaining quantities for internal layers:
            cosha[i], sinha[i] = sincos(im * layers[i].width * γ)
            sinha[i] = -sinha[i]
            tanha[i] = sinha[i] / cosha[i]
            Zleft[i] = Z0[i] *
                (Zleft[i-1] + Z0[i]*tanha[i]) / (Z0[i] + Zleft[i-1] * tanha[i])  # Eq. (6.12b).
        end
        # Loop over regions in backwards order to get Zright:
        Zright[s] = zero(eltype(Zright)) # short-circuited
        #Zright[s-1] = Z0[s] * tanha[s] # Eq. (6.23a)
        for i in s-1:-1:1
            Zright[i] = Z0[i+1] * (Zright[i+1] + Z0[i+1]*tanha[i+1]) / 
                (Z0[i+1] + Zright[i+1] * tanha[i+1])  # Eq. (6.23b)
        end

        # Compute voltage and current at z = z_1 using Eq. (6.22):
        current = 2 / (1 + Zright[1] / Z0[1])
        voltage = current * Zright[1]
        # Compute reflection term using Eq. (6.27):
        gsm.s11[q1,q1] = 1 - current
      
        # Step voltage and current up to junction s using Eq. (6.24):
        for i in 2:s
            current = current * cosha[i] - voltage * sinha[i] / Z0[i]
            voltage = current * Zright[i]
        end
        iincs[q1,1] = current  # Save incident current.

      
        # Now compute TLGF:
        current =  1 / Zleft[s] # Eq. (6.38a)
        voltage = -one(current)     # Eq. (6.38a)
        # Step down to z_1 using (6.32):
        for i in s-1:-1:1
            current = current * cosha[i+1] + voltage * sinha[i+1] / Z0[i+1] 
            voltage = -current * Zleft[i]
        end
        tlgfiv[q1,1] = current # Save quantity used in Eq. (6.37a)

    end  # Done with Region 1 modes
    
    #  Loop over the region N source modes:
    for q2 in 1:n2
        Z0[N] = 1 / layers[N].Y[q2]  # Obtain region N modal impedance
        Zright[N-1] = Z0[N]  # Semi-infinite line
        # Compute mag squared of transverse wavenumber (rad^2/m^2):
        β² = layers[N].β[q2] ⋅ layers[N].β[q2] 
        for i in N-1:-1:s+1 # Loop over each region
            # Obtain wavenumber (Region i) squared (rad^2/m^2):
            k² = k0sq * layers[i].μᵣ * layers[i].ϵᵣ
            γ = mysqrt(β² - k²) # Atten. constant (neper/m).
            if layers[N].P[q2] == TE
                Z0[i] =  (im * (k0*η₀) * layers[i].μᵣ) / γ # TE
            else
                Z0[i] = γ / (im * k0/η₀ * layers[i].ϵᵣ)
            end
            # Compute remaining quantities for internal layers:
            cosha[i], sinha[i] = sincos(im * layers[i].width * γ)
            sinha[i] = -sinha[i]
            tanha[i] = sinha[i] / cosha[i]
            Zright[i-1] = Z0[i] * (Zright[i] + Z0[i]*tanha[i]) / 
                (Z0(i) + Zright(i) * tanha(i))  # Eq. (6.23b)
        end
        # Loop over regions in forward order to get Zleft:
        Zleft[s] = zero(eltype(Zleft)) # Short-circuited
        #Zleft(s+1) = Z0(s+1) * tanha(s+1) # Eq. (6.31a)
        for i in s+1:N-1
            Zleft[i] = Z0[i] * (Zleft[i-1] + Z0[i]*tanha[i]) / 
                 (Z0[i] + Zleft[i-1] * tanha[i])  # Eq. (6.31b)
        end

        # Compute voltage and current at z = z_{N-1} using Eq. (6.30):
        current = -2  / (1 + Zleft[N-1] / Z0[N])
        voltage = -current * Zleft[N-1]
        # Compute reflection term using Eq. (6.35):
        gsm.s22[q2,q2] = 1 + current
      
        # Step voltage and current down to junction s using Eq. (6.32):
        for i in N-2:-1:s
            current = current * cosha[i+1] + voltage * sinha[i+1] / Z0[i+1]
            voltage = -current * Zleft[i]
        end
        iincs[q2,2] = current  # Save incident voltage at FSS plane

        # Now compute TLGF's:
        current = 1 / Zright[s] # Eq. (6.38b)
        voltage = one(current) # Eq. (6.38b)
        # Step up to z_{N-1} using (6.24):
        for i in s+1:N-1
            current = current * cosha[i] - voltage * sinha[i] / Z0[i]
            voltage = current * zright[i]
        end
        tlgfiv[q2,2] = current # Save quantity used in (6.18b)
    end # Done with Region N modes
    return (gsm, tlgfiv, iincs)
end



"""
    find_mode_index(p::TEorTM, m::Int, n::Int, layer::Layer)
    
Find the index of the mode (p,m,n).

## Arguments

- `p`, `m`, `n`: Mode indices. 
- `layer`: An instance of Layer with modes initialized.

## Return Value

The positive integer index, if found.  Otherwise, zero.
"""
function find_mode_index(p::TEorTM, m::Int, n::Int, layer::Layer)
    for i in 1:length(layer.P)
        if layer.P[i] == p  && layer.M[i] == m && layer.N[i] == n
            return i
        end
    end
    return 0
end




"""
    gsm_slab_interface(L1::Layer, L2::Layer, k0) -> gsm::GSM

Calculate the GSM (generalized scattering matrix) of a dielectric junction.

## Arguments

- `L1`, `L2`:  Contains the layer parameters layers 1 and 2, resp.
- `k0`: Free-space wavenumber in rad/m.

## Return Values

- `gsm`: A variable containing the generalized scattering matrix computed in accordance
    with Equations (3.23) and (3.24) of the theory documentation.

"""
function gsm_slab_interface(L1::Layer, L2::Layer, k0)::GSM

      n1 = length(L1.P)  # Number of modes in Region 1.
      n2 = length(L2.P)  # Number of modes in Region 2.
      gsm = GSM(n1, n2)
      #  Loop over the region 1 source modes:
      for i1 in 1:n1
          Y1 = L1.Y[i1]  # Obtain source modal admittance.
          # We now require the modal admittance for the corresponding
          # mode in Region 2. This, however, is not guaranteed to exist 
          # since we don't necessarily define the same sets of modes in all 
          # regions.  Therefore, we search for it on the other side,
          # but if it doesn't exist, we have to compute the modal 
          # admittance on the fly.
      
          # Locate identical mode in region 2:
          i2 = find_mode_index(L1.P[i1], L1.M[i1], L1.N[i1], L2)
          if i2 == 0
              # Obtain wavenumber (Region 2) squared
              k2sq = L2.β[1] ⋅ L2.β[1] - L2.γ[1]^2
              β² = L1.β[i1] ⋅ L1.β[i1] # for desired mode
              γ = mysqrt(β² - k2sq)
              if L1.P[i1] == TE 
                  Y2 = γ / (im * (k0*η₀) * L2.μᵣ)
              else
                  Y2 = im * (k0/η₀) * L2.ϵᵣ / γ
              end
          else
              Y2 = L2.Y[i2]
          end

          # Compute reflection and transmission coefficients:
          R = (Y1 - Y2) / (Y2 + Y1) # Eq. (3.24a)
          gsm.s11[i1,i1] = R
          if i2 != 0
              T = 2 * mysqrt(Y2) *  mysqrt(Y1) / (Y2 + Y1) # Eq. (3.24b)
              gsm.s21[i2,i1] = T
          end
      end

      #  Loop over the region 2 source modes:
      for i2 in 1:n2
          Y2 = L2.Y[i2]  # Obtain source modal admittance.
          # We now require the modal admittance for the corresponding
          # mode in Region 1. This, however, is not guaranteed to exist 
          # since we don't necessarily define the same sets of modes in all 
          # regions.  Therefore, we search for it on the other side,
          # but if it doesn't exist, we have to compute the modal 
          # admittance on the fly.
      
          # Locate identical mode in region 1:
          i1 = find_mode_index(L2.P[i2], L2.m[i2], L2.n[i2], L1)

          if i1 == 0
              # Obtain wavenumber (Region 1) squared
              k1sq = L1.β[1] ⋅ L1.β[1] - L1.γ[1]^2
              β² = L2.β[i2] ⋅ L2.β[i2] # for desired mode
              γ = mysqrt(β² - k1sq)
              if L2.P[i2] == TE
                  Y1 = γ / (im * (k0*η₀) * L1.μᵣ)
              else
                  Y1 = im * (k0/η₀) * L1.ϵᵣ / γ
              end
          else
              Y1 = L1.Y[i1]
          end
          # Compute reflection and transmission coefficients:
          R = (Y1 - Y2) / (Y2 + Y1) # Eq. (3.24a)
          gsm.s22[i2,i2] = -R
          if i1 != 0
              T = 2 * mysqrt(Y2) * mysqrt(Y1) / (Y2 + Y1) # Eq. (3.24b)
              gsm.s12[i1,i2] = T
          end
      end
      return gsm
  end # function

"""
    initialize_gsm_file(fname::AbstractString, layer1::Layer, layer2::Layer, sh1::Sheet, sh2::Sheet) 

Initialize a generalized scattering matrix (GSM) file, by appending the GSM-specific 
header info, including the unit cell data.  The file is stored in `JLD2` format.


## Arguments

- `fname`: A string containing the file name to open.  The file will be closed after writing.
- `layer1`, `layer2`:  These provide the mode indices for the input and output regions.
- `sh1`, `sh2`: `SHEET` objects from which the unit 
         cell info for the input and output regions will be extracted.  These two arguments
    are optional, but if either is present the other must be also.  If they are absent, then
    a square unit cell of dimension 1 meter will be used for both input and output regions.
"""
function initialize_gsm_file(fname::AbstractString, layer1::Layer, layer2::Layer, sh1::Sheet, sh2::Sheet)
    jldopen(fname, "w") do fid
        fid["FileType"] = "GSMFile"
        fid["FileTypeVersion"] = v"0.0.1"
        # Unit cell information:
        fid["s1in"] = sh1.s₁
        fid["s2in"] = sh1.s₂
        fid["unitsin"] = sh1.units
        fid["s1out"] = sh2.s₁
        fid["s2out"] = sh2.s₂
        fid["unitsout"] = sh2.units
        
        # Modelists:
        fid["ModelistIn"] = (layer1.P, layer1.M, layer1.N)
        fid["ModelistOut"] = (layer2.P, layer2.M, layer2.N)
    end
end # function

function initialize_gsm_file(fname::AbstractString, layer1::Layer, layer2::Layer)
    sh = RWGSheet()
    sh.s₁ = [1.0, 0.0]
    sh.s₂ = [0.0, 1.0]
    sh.units = u"m"
    initialize_gsm_file(layer1, layer2, fname, sh, sh)
    return
end


"""
    append_gsm_data(fname::AbstractString, gname::String, gsm::GSM, l1::Layer, l2::Layer, case::Dict{String,Float64})

Append a block of GSM data to a GSM file for a particular frequency and pair of scan parameters.  It is 
assumed that the file has already been initialized by a previous call to `initialize_gsm_file`.

## Arguments

- `fname`: The name of the GSM file to be appended to.  The file should already have been initialized
via a call to `initialize_gsm_file`.
- `gname`: The unique `JLD2` group name to be used in the file for grouping the data associated with this frequency/scan case.
- `gsm`:  The GSM data to be written to the file.
- `l1`, `l2`: `Layer` data for the input and output regions containing the permeability and permittivity information
to be written to the file.
- `case`: A dictionary containing three key-value pairs which together fully define the frequency/scan case. 
The keys are "FGHz" and either ("θ" and "ϕ") or ("ψ₁" and "ψ₂").
"""
function append_gsm_data(fname::AbstractString, gname::String, gsm::GSM, l1::Layer, l2::Layer, case::Dict{String,Float64})
    jldopen(fname, "a") do fid
        group = JLD2.Group(fid, gname)
        group["case"] = case
        group["epsmu1"] = (l1.ϵᵣ, l1.μᵣ)
        group["epsmu2"] = (l2.ϵᵣ, l2.μᵣ)
        group["gsm"] = gsm
    end
    return    
end

function read_gsm_file(fname::AbstractString)
    dat = load(fname)
    (haskey(dat, "FileType") && dat["FileType"] == "GSMFile") || error("$fname is not a GSM file")
    dat
end
    
"""
    translate_gsm!(g::GSM, dx::Real, dy::Real, layer1::Layer, layer2::Layer)

Convert `GSM` `g` into one for an identical structure that is translated by 
`dx` meters in x and `dy` meters in y. `layer1` and `layer2` define the modes 
on each side of the GSM.

## Arguments:

- `g`: On input contains the GSM for a device with zero x and y translation. After
   return, `g` will be modified by applying appropriate phase shifts due to translation,
   according the to the formula presented in John C. Vardaxoglou, "Frequency Selective Surfaces", 
   Wiley, 1997, Eq. (7.57).
- `dx`, `dy`:   The translations in the x and y directions in meters.
- `layer1`, `layer2`:  The β components of these variables are used with 
    `dx` and `dy` to determine the effect of the translation.
"""
function translate_gsm!(g, dx, dy, layer1, layer2)
    
    COMPLEX(WP)                          :: t1(size(layer1.p)), t2(size(layer2.p))
    dxy = SA[dx, dy]

    # Check that layers have correct number of modes:
    size(g.s11,1) == length(layer1.P) || error("N1 mismatch")
    size(g.s21,1) == length(layer2.P) || error("N2 mismatch")

    # Compute phase shift constants:
    t1 = [cis(β ⋅ dxy) for β in layer1.β]
    t2 = [cis(β ⋅ dxy) for β in layer2.β]

    # Fix up S11:
    for n in 1:length(layer1.P), m in 1:length(layer1.P)
        g.s11[m,n] *= t1[n] / t1[m]
    end

    # Fix up S12:
    for n in 1:length(layer2.P), m in 1:length(layer1.P)
        g.s12[m,n] *= t2[n] / t1[m]
    end

    # Fix up S21:
    for n in 1:length(layer1.P), m in 1:length(layer2.P)
        g.s21[m,n] *= t1[n] / t2[m]
    end

    # Fix up S22:
    for m in 1:length(layer2.P), n in 1:length(layer2.P)
        g.s22[m,n] *= t2[n] / t2[m]
    end

    return
end


"""
    choose_gblocks(strata::Vector{Union{Layer,Sheet}}, k0min) -> gbl::Vector{Gblock}

Set up gbl, the array of GBLOCKs the defines the basic GSM building 
blocks for the FSS structure.

## Input arguments:

- `strata`: A vector whose elements are either `Layer` or `Sheet` objects. `strata` defines the planar
geometry being analyzed.
- `k0min`:    The minimum free-space wavenumber (rad/m) for all analysis cases.

## Return value:

- `gbl`:  A vector of `Gblock` elements defining the GBLOCKS used to analyze the complete FSS structure.
"""
function choose_gblocks(strata::Vector{Union{Layer,Sheet}}, k0min)::Vector{Gblock}
    islayer = map(t -> t isa Layer, strata)
    layers = @view strata[islayer]
    nl = length(layers)
    nj = nl-1 # number of dielectric junctions
    issheet = map(x -> x isa Sheet, strata)
    sheets = @view strata[issheet]
    ns = length(sheets)
    if issheet[1] || issheet[end]
        error("First and last strata must be dielectric layers")
    end
    if ns > 1 && any(diff(findall(issheet)) .== 1)
        error("Adjacent sheets must be separated by one or more dielectric layers")
    end

    # Pure radome case is easy:
    if ns == 0 || all(t.style == "NULL" for t in sheets)
        gbl = [Gblock(i:i, 0) for i = 1:nj]
        return gbl
    end

    # Compute min. electrical length of each layer:
    elength = [k0min * l.width * sqrt(real(l.ϵᵣ * l.μᵣ)) for l in layers]

    sint = cumsum(islayer)[issheet] # sint[k] contains dielectric interface number of k'th sheet 

    junc = zeros(Int, nj)
    junc[sint] = 1:ns #  junc[i] is the sheet number present at interface i, or 0 if no sheet is there

    upa = find_unique_periods(junc, sheets)
    
    # Initialize junction ownerships to self:
    owner = collect(1:nj) # owner[i] is the interface number to which interface i belongs.
      
    # Assign junction ownership for thin layers adjacent and prior to location of first sheet:
    if sheets[1].style ≠ "NULL" 
        for j in sint[1]:-1:2 # layer index
            elength[j] > min_elength && break
            owner[j-1] = sint[1]
        end
    end
    # Assign junction ownership for thin layers adjacent to and after location of last sheet:
    if sheets[end].style ≠ "NULL" 
        for j in 1+sint[ns]:nj # layer index
            elength[j] > min_elength && break
            owner[j] = sint[end]
        end
    end
      
    # Loop over each adjacent pair of junctions stored in sint:
    for i in 1:ns-1
        (j1, j2) = sint[i:i+1]   # Interface locations
        ## to do:  check here if j2-j1==1, and periods are different
        Nwide = count(elength[(1+j1):j2] .> min_elength) # Numb. of wide enough intervening layers
        if Nwide > 0  # step into the region from each end and assign ownership
            for j in j2:-1:j1+1
                (elength[j] > min_elength || sheets[junc[j2]].style == "NULL") && break
                owner[j-1] = j2
            end
            for j in 1+j1:j2
                (elength[j] > min_elength || sheets[junc[j1]].style == "NULL") && break
                owner[j] = j1
            end
        else # here Nwide==0, so none of the layers are really wide enough
            @warn """
                    Layers are too thin between sheets $(junc[j1]) and $(junc[j2]) at interfaces $(j1) and $(j2).
                             GSM calculation may not be accurate.
                             Continuing with execution anyway."""
            
            # Assign neighboring junctions for thin layers adjacent to j1 and j2 to j1 and j2
            #
            if sheets[junc[j1]].style == "NULL" && sheets[junc[j2]].style ≠ "NULL"
                owner[j1+1:j2] .= j2
                continue
            elseif sheets[junc[j1]].style ≠ "NULL" && sheets[junc[j2]].style == "NULL"
                owner[j1:j2-1] .= j1
                continue
            elseif sheets[junc[j1]].style == "NULL" && sheets[junc[j2]].style == "NULL"
                continue
            end
            #
            # At this point, we've determined that both sheets are non-null, and closely separated
            #
            jmid = (j1 + j2 + 1) ÷ 2
            if 2jmid > j1+j2 && upa[j1] ≠ upa[j2] # odd number of thin layers between sheets of different periodicities
                error("Inconsistent periodicities for sheets $(junc[j1]) and $(junc[j2]). Try splitting layer $(jmid)")
            end
            # We have an even number of thin layers between.  We will not include the center two layers in the gblocks
            # so that we can match up to the modes in the gblocks on either side.
            owner[j1:jmid-1] .= j1
            owner[jmid+1:j2] .= j2
        end
    end
      
    # Count number of GBLOCKS needed:
    ngbl = 1
    for j in 2:nj
        owner[j] ≠ owner[j-1] && (ngbl += 1)
    end
    # and assign them:
    gbl = [Gblock(1:1, 0) for j in 1:nglb]
    jg = 1 # Initialize Gblock counter
    for j in 1:nj  # Loop over each interface
        if owner[j] == owner[first(gbl[jg].rng)]
            gbl[jg].rng = first(gbl[jg].rng):j # Increment the end pointer.
        else
            jg += 1
            gbl[jg].range = j:j
            gbl[jg].j = 0
        end
        junc[j] ≠ 0 && (gbl[jg].j = j)
    end

    return gbl
end

end # module
