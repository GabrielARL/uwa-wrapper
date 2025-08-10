function [range, depths, TL] = compute_transmisson_loss(args)
% COMPUTE_TRANSMISSON_LOSS
% Always write input_structs.mat from args, run Julia, and return 2D TL grid.
% Outputs:
%   range  : 1×Nr (m)
%   depths : 1×Nz (m, positive down)
%   TL     : Nz×Nr (dB), as saved by transmission_loss_kraken.jl

% --- Write input_structs.mat from args ---
svp_in      = args.svp;
env_in      = args.env;
opt_in      = args.opt;
output_mode = struct( ...
    'rmin',       args.out.rmin, ...
    'rmax',       args.out.rmax, ...
    'range_step', args.out.dr, ...
    'depth_step', args.out.dz, ...   % ensure args.out.dz is set in run_sim_TL.m
    'zm',         args.out.zm );

save('input_structs.mat', 'svp_in', 'env_in', 'opt_in', 'output_mode', '-v7');

% --- Call Julia wrapper ---
cmd = 'env -i HOME=$HOME PATH=$PATH LC_ALL=C julia --project=. transmission_loss_kraken.jl';
status = system(cmd);
if status ~= 0
    error('❌ Julia transmission_loss_kraken.jl failed (status=%d).', status);
end

% --- Load outputs ---
D       = load('pekeris_kraken_output.mat');
range   = D.range;      % 1×Nr
depths  = D.depths;     % 1×Nz
TL      = D.tl_field;   % Nz×Nr (dB)
end
