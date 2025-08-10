% run_sim_TL.m — Example driver for the Kraken TL wrapper
clear; clc;

args = struct();

% --- SVP ---
args.svp.c0 = 1500;  % m/s
args.svp.wssp = [0 1500]; % minimal SSP table, depth (m) vs speed (m/s)

% --- Environment & seabed ---
args.env.bathymetry = 50;          % m (scalar)
args.env.water_density = 1000;       % kg/m^3
args.env.seabed.c = 2000;            % m/s
args.env.seabed.rho = 2000;          % kg/m^3
args.env.seabed.alpha = 0.5;         % attenuation (dummy)

% --- Seabed layering (dummy support) ---
% Format: [layer_index, thickness(m), c(m/s), rho(kg/m^3)]
args.env.seabed_layer = [
    1,  40, 1520, 1540;
    2,  50, 1730, 1730
];

% --- Operating point ---
args.opt.fcw = 300;   % Hz
args.opt.ztx = 10;   % m (positive down)

% --- Output grid ---
args.out.rmin = 0;     % m
args.out.rmax = 1000000;  % m
args.out.dr   = 1000;    % m
args.out.dz   = 1;     % m   % vertical resolution
args.out.zm   = 25;    % m (positive depth for TL line output)

[range, depths, TL] = compute_transmisson_loss(args);

% Convert to dB (20 log10)
TL_db = 20*log10(abs(TL));

figure('Color','w');
imagesc(range, depths, TL_db);
set(gca,'YDir','normal');
xlabel('Range (m)');
ylabel('Depth (m)');
title(sprintf('Pekeris (Kraken) — TL (f = %g Hz)', args.opt.fcw));
colormap(turbo);
hcb = colorbar;
ylabel(hcb, 'TL (dB)');
