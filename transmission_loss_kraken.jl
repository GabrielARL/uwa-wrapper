# transmission_loss_kraken.jl — 2D TL grid with depth-dependent SSP (SampledField-aware)
using MAT
using UnderwaterAcoustics
using AcousticsToolbox
using Printf
using Statistics

# ---------- Helpers ----------

"""
    build_soundspeed(svp_in) -> Union{Float64, SampledField}

Accepts:
  - svp_in["wssp"] :: Nx2 [depth(m, +down), c(m/s)]
      -> SampledField with :cubic interpolation if depths are evenly spaced,
         otherwise :linear.
  - svp_in["c0"]   :: scalar m/s
Depths are positive downward; UA uses z≤0 downward, so we map d -> z = -d.
"""
function build_soundspeed(svp_in)
    if haskey(svp_in, "wssp")
        w = svp_in["wssp"]
        @assert size(w,2) == 2 "wssp must be Nx2 [depth(m), c(m/s)]"
        d = Float64.(w[:,1])                 # positive depths
        c = Float64.(w[:,2])                 # speeds
        p = sortperm(d); d = d[p]; c = c[p]  # sort by depth
        zvec = -d                            # z ≤ 0 is down

        # If evenly spaced, allow :cubic (requires AbstractRange)
        uniform = length(zvec) ≥ 2 ? all(abs.(diff(zvec) .- (zvec[2]-zvec[1])) .< 1e-9) : false
        if uniform
            Δ  = zvec[2] - zvec[1]
            zr = range(zvec[1]; step=Δ, length=length(zvec))  # AbstractRange
            return SampledField(c; z=zr, interp=:cubic)
        else
            return SampledField(c; z=zvec, interp=:linear)
        end
    elseif haskey(svp_in, "c0")
        return Float64(svp_in["c0"])
    else
        return 1500.0
    end
end

"""
    normalize_bathymetry(bathy_any) -> Float64

Scalar bathymetry passes; Nx2 [range, depth] reduced to mean depth (dummy for Kraken input).
"""
function normalize_bathymetry(bathy_any)::Float64
    if isa(bathy_any, Number)
        return Float64(bathy_any)
    else
        A = Float64.(bathy_any)
        return mean(A[:,2])
    end
end

"""
    make_bathy_func(bathy_any) -> (r::Real) -> depth(m)

Scalar -> constant; Nx2 -> piecewise linear, clamped at ends.
"""
function make_bathy_func(bathy_any)::Function
    if isa(bathy_any, Number)
        H = Float64(bathy_any)
        return r -> H
    else
        A = Float64.(bathy_any)
        rs, hs = A[:,1], A[:,2]
        if !issorted(rs)
            p = sortperm(rs); rs = rs[p]; hs = hs[p]
        end
        function depth_at(r::Real)
            x = Float64(r)
            if x <= rs[1]; return hs[1]; end
            if x >= rs[end]; return hs[end]; end
            i = searchsortedlast(rs, x); i = min(max(i,1), length(rs)-1)
            t = (x - rs[i])/(rs[i+1]-rs[i])
            return (1-t)*hs[i] + t*hs[i+1]
        end
        return depth_at
    end
end

"""
    make_seabed(sb) -> FluidBoundary
Requires sb["c"], sb["rho"]; optional sb["alpha"].
"""
function soundspeed_from_svp(svp_in)
    if haskey(svp_in, "wssp")
        w = svp_in["wssp"]
        @assert size(w,2) == 2 "wssp must be Nx2 [depth, c]"
        d = Float64.(w[:,1])    # depths, +down
        c = Float64.(w[:,2])    # speeds
        p = sortperm(d); d = d[p]; c = c[p]
        z = -d                  # convert to model z (≤0 is down)

        # evenly spaced? -> AbstractRange for :cubic
        uniform = length(z) ≥ 2 ? all(abs.(diff(z) .- (z[2]-z[1])) .< 1e-9) : false
        if uniform
            Δz = z[2]-z[1]
            zrange = range(z[1]; step=Δz, length=length(z))  # AbstractRange
            return SampledField(c; z=zrange, interp=:cubic)
        else
            return SampledField(c; z=z, interp=:linear)
        end
    elseif haskey(svp_in, "c0")
        return Float64(svp_in["c0"])
    else
        return 1500.0
    end
end

function make_seabed(sb)::FluidBoundary
    c_sb   = Float64(sb["c"])
    rho_sb = Float64(sb["rho"])
    if haskey(sb, "alpha")
        try
            return FluidBoundary(c_sb, rho_sb, Float64(sb["alpha"]))
        catch
            return FluidBoundary(c_sb, rho_sb)
        end
    else
        return FluidBoundary(c_sb, rho_sb)
    end
end

# ---------- Load inputs ----------
mat = matopen("input_structs.mat")
svp_in      = read(mat, "svp_in")
env_in      = read(mat, "env_in")
opt_in      = read(mat, "opt_in")
output_mode = read(mat, "output_mode")
close(mat)

# ---------- Scalars, fields, and grids ----------
ρw   = Float64(env_in["water_density"])
sb   = make_seabed(env_in["seabed"])
cw   = build_soundspeed(svp_in)                  # Float64 or SampledField

fc   = Float64(opt_in["fcw"])
z_tx = haskey(opt_in, "ztx") ? -abs(Float64(opt_in["ztx"])) : -10.0

rmin = Float64(output_mode["rmin"])
rmax = Float64(output_mode["rmax"])
dr   = Float64(output_mode["range_step"])
dz   = Float64(get(output_mode, "depth_step", 1.0))
zpos = Float64(output_mode["zm"])                # kept for 1D compatibility
z_rx = -abs(zpos)

# StepRange for x (range)
xline = rmin:dr:rmax

# Bathymetry function and checks
bathy_fn = make_bathy_func(env_in["bathymetry"])
bathys   = map(bathy_fn, collect(xline))
rx_depth = abs(z_rx)
tx_depth = abs(z_tx)

invalid_rx = findall(d -> rx_depth > d + 1e-9, bathys)
if !isempty(invalid_rx)
    bad_rs = join(round.(collect(xline)[invalid_rx[1:min(end,10)]]; digits=3), ", ")
    error("Receiver depth $(rx_depth)m exceeds bathymetry at $(length(invalid_rx)) ranges. First: [$bad_rs]")
end
H_tx = bathy_fn(0.0)
if tx_depth > H_tx + 1e-9
    error("Transmitter depth $(tx_depth)m exceeds bathymetry $(H_tx)m at range 0 m")
end

# Depth grid as StepRange (negative z), and positive-down copy for MATLAB
H_min  = minimum(bathys)
zgrid  = (-H_min):dz:0.0                   # StepRange (z ≤ 0)
depths = collect(0.0:dz:H_min)             # positive-down for MATLAB

# ---------- Environment & model ----------
# Kraken often expects numeric cw in its .env writer; be conservative:
c_eff = cw isa Number ? cw :
        (haskey(svp_in,"wssp") ? mean(Float64.(svp_in["wssp"][:,2])) : 1500.0)

H_mean = normalize_bathymetry(env_in["bathymetry"])     # RD->mean for Kraken
cw = soundspeed_from_svp(svp_in)  # Number or SampledField
env = UnderwaterEnvironment(
    bathymetry = H_mean,
    soundspeed = cw,                                  # scalar for Kraken stability
    density    = ρw,
    seabed     = sb,
)

maxc = haskey(svp_in,"wssp") ? maximum(Float64.(svp_in["wssp"][:,2])) : c_eff
pm = Kraken(env; chigh = 1.1*max(1500.0, maxc))

tx = AcousticSource(0.0, z_tx, fc)

# ---------- 2D TL (dB) ----------
rxs   = AcousticReceiverGrid2D(xline, zgrid)   # StepRanges on both axes
TL_db = transmission_loss(pm, tx, rxs)         # size: length(zgrid) × length(xline)

# ---------- Save outputs ----------
mat_out = matopen("pekeris_kraken_output.mat", "w")
write(mat_out, "range", collect(xline))
write(mat_out, "depths", depths)
write(mat_out, "tl_field", TL_db)              # already in dB from Kraken wrapper
close(mat_out)
