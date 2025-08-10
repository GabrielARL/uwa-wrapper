# transmission_loss_kraken.jl — 2D TL grid (saves depths & tl_field)
using MAT
using UnderwaterAcoustics
using AcousticsToolbox
using Printf
using Statistics

# ---- Helpers ----
function effective_c(svp_in)::Float64
    if haskey(svp_in, "c0")
        return Float64(svp_in["c0"])
    elseif haskey(svp_in, "wssp")
        return mean(Float64.(svp_in["wssp"][:, 2]))
    else
        return 1500.0
    end
end

function normalize_bathymetry(bathy_any)::Float64
    if isa(bathy_any, Number)
        return Float64(bathy_any)
    else
        A = Float64.(bathy_any)
        return mean(A[:,2])
    end
end

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

# ---- Load inputs ----
mat = matopen("input_structs.mat")
svp_in      = read(mat, "svp_in")
env_in      = read(mat, "env_in")
opt_in      = read(mat, "opt_in")
output_mode = read(mat, "output_mode")
close(mat)

# ---- Scalars & functions ----
H_mean = normalize_bathymetry(env_in["bathymetry"])
ρw     = Float64(env_in["water_density"])
sb     = make_seabed(env_in["seabed"])

fc   = Float64(opt_in["fcw"])
z_tx = haskey(opt_in, "ztx") ? -abs(Float64(opt_in["ztx"])) : -10.0

rmin = Float64(output_mode["rmin"])
rmax = Float64(output_mode["rmax"])
dr   = Float64(output_mode["range_step"])
dz   = Float64(get(output_mode, "depth_step", 1.0))
zpos = Float64(output_mode["zm"])            # keep for compatibility
z_rx = -abs(zpos)

# Range grid
ranges = collect(rmin:dr:rmax)

# Bathymetry checks (RX at all ranges; TX at r=0)
bathy_fn = make_bathy_func(env_in["bathymetry"])
rx_depth = abs(z_rx); tx_depth = abs(z_tx)
bathys   = map(bathy_fn, ranges)

invalid_rx = findall(d -> rx_depth > d + 1e-9, bathys)
if !isempty(invalid_rx)
    bad_rs = join(round.(ranges[invalid_rx[1:min(end,10)]]; digits=3), ", ")
    error("Receiver depth $(rx_depth)m exceeds bathymetry at $(length(invalid_rx)) ranges. First: [$bad_rs]")
end
H_tx = bathy_fn(0.0)
if tx_depth > H_tx + 1e-9
    error("Transmitter depth $(tx_depth)m exceeds bathymetry $(H_tx)m at range 0 m")
end

# Depth grid (0..min depth along track)
H_min = minimum(bathys)
Nz = max(1, Int(floor(H_min/dz)) + 1)
depths = range(0.0, step=dz, length=Nz) |> collect  # positive down

# ---- Environment & model ----
c_eff = effective_c(svp_in)
env = UnderwaterEnvironment(
    bathymetry = H_mean,      # single representative depth for Kraken input
    soundspeed = c_eff,       # scalar (Kraken expects numeric)
    density    = ρw,
    seabed     = sb,
)

maxc = haskey(svp_in,"wssp") ? maximum(Float64.(svp_in["wssp"][:,2])) : c_eff
pm = Kraken(env; chigh = 1.1*max(1500.0, maxc))

tx = AcousticSource(0.0, z_tx, fc)

# --- Build grids as StepRanges (NOT vectors) ---
xline = rmin:dr:rmax                  # StepRangeLen, good
H_min = minimum(bathys)
zgrid = (-H_min):dz:0.0               # negative (down) up to 0, also StepRangeLen

# for saving to MATLAB (positive depth vector for imagesc)
depths = collect(0.0:dz:H_min)

# 2D TL
rxs    = AcousticReceiverGrid2D(xline, zgrid)
TL_db  = transmission_loss(pm, tx, rxs)   # size: length(zgrid) × length(xline)

# single-depth line (keep using the StepRange for xline)
rxs_ln   = AcousticReceiverGrid2D(xline, z_rx)
TL_line  = transmission_loss(pm, tx, rxs_ln) |> vec

# save
mat_out = matopen("pekeris_kraken_output.mat", "w")
write(mat_out, "range", collect(xline))   # okay to collect for MATLAB I/O
write(mat_out, "depths", depths)
write(mat_out, "tl_field", TL_db)
write(mat_out, "tl_line", TL_line)
close(mat_out)

