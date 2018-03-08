% This function provides an estimated depth scale for radar images based 
% off of depth-density relations from nearby firn cores

function [radar, core_composite] = radar_depth(radar_file, cores)

% Firn core data must be provided in the proper format (i.e. as outputted
% by "import_cores.m"

% Import, clean, and reformat radar data
radar = radar_clean(radar_file);

% Import and clean radar data, and estimate depth-density profile in the 
% radar based on spatial weighting of nearby cores
[core_composite] = rho_spW(radar, cores);

% Needed constants
rho_ice = 0.917;    % g/cm^3
rho_0 = mean(core_composite.rho(1:10));
emiss_ice = 3.17;

% g = @(a, b, c, x) a*(x).^b + c;
% radar.rho_coeff = nlinfit(core_composite.depth, core_composite.rho, g,...
%     [1; 0.5; core_composite.rho(1)]);
% depth_mod = (0:0.02:150)';
% rho_mod = radar.rho_coeff(1)*depth_mod.^radar.rho_coeff(2) + radar.rho_coeff(3);

% % Fit depth-density model (exponential function) to synthetic core data
% g = @(m, x) rho_ice - (rho_ice-rho_0)*exp(m*x);
% radar.rho_coeff = nlinfit(core_composite.depth, core_composite.rho, g, 0);
% % Calculate modeled density-depth data based on fitted exp function
% depth_mod = (0:0.02:150)';
% rho_mod = rho_ice - (rho_ice-rho_0)*exp(radar.rho_coeff*depth_mod);



% % Fit depth-density model (exponential function) to synthetic core data
% % g = @(s, b, m, c, x) rho_ice + s*b.^(m*x+c);
% g = @(m, x) rho_ice - (rho_ice-rho_0)*exp(m*x);
% radar.rho_fit = fit(core_composite.depth, core_composite.rho, g,...
%     'StartPoint', 0);
% % Calculate modeled density-depth data based on fitted exp function
% depth_mod = (0:0.02:150)';
% rho_coeff = coeffvalues(radar.rho_fit);
% rho_mod = rho_ice - (rho_ice-rho_0)*exp(rho_coeff*depth_mod);

g = @(a, b, c, x) a*(x).^b + c;
radar.rho_fit = fit(core_composite.depth, core_composite.rho, g,...
    'StartPoint', [1 0.5 core_composite.rho(1)]);
depth_mod = (0:0.02:150)';
rho_coeff = coeffvalues(radar.rho_fit);
rho_mod = rho_coeff(1)*depth_mod.^rho_coeff(2) + rho_coeff(3);
% 
% Diagnostic figures for density modeling
% figure
% hold on
% plot(core_composite.rho, core_composite.depth)
% plot(rho_mod, depth_mod)
% ylim([0 core_composite.depth(end)])
% ylabel('Depth (m)')
% xlabel('Density (g/cm^3)')
% set(gca, 'Ydir', 'reverse')
% legend('Composite firn core', 'Density model')
% hold off

% Calculate the real part of the relative permittivity of firn
e_mod = (1 + 0.845*rho_mod).^2;  % Kovacs
% e_mod2 = ((rho_mod/rho_ice)*(emiss_ice^(1/3)-1) + 1).^3; % Looyenga, 1965

% Calculate radar propagation speed with depth from emissivity
c0 = 2.9979E8;  % speed of light in a vacuum (m/s)
cZ_mod = c0./sqrt(e_mod);   % Speed of light with depth in firn (m/s)

% Calculate modeled one way travel time (seconds)
time_disc = mean(diff(depth_mod))./cZ_mod;
time_mod = cumsum([0; time_disc(1:end-1)]);
% time_mod = cumsum([0; diff(depth_mod)]./cZ_mod);

% Calculate cumulative two-way travel time for radar data
radar.TWTT = (0:radar.time_trace(1):radar.time_trace(1)*...
        (size(radar.data_out, 1)-1))';

% Generate array of cumulative two-way travel time for radar data based on
% recorded TWT time (for 2010 data) or recorded one-way travel time (for
% 2011 data)
% if contains(radar_file, '2010') == true
%     radar.TWTT = (cumsum(2*(0:radar.time_trace(1):radar.time_trace(1)*...
%         (size(radar.data_out, 1)-1))))';
% else
%     radar.TWTT = (cumsum(2*(0:radar.time_trace(1):radar.time_trace(1)*...
%         (size(radar.data_out, 1)-1))))';
% end

% Convert recorded TWT time to depth by interpolating measured values to
% the modeled time-depth relationship
radar.depth = interp1(time_mod, depth_mod, 0.5*radar.TWTT);



% % Calculate radar depth based off of raw (no model) density data in
% % synthetic core (for error calculations)
% e_data = (1 + 0.845*core_synthetic.rho).^2;  % Kovacs
% % e_data2 = ((core_synthetic.rho/rho_ice)*(emiss_ice^(1/3)-1) + 1).^3; % Looyenga, 1965
% cZ_data = c0./sqrt(e_data);
% time_data = cumsum(core_synthetic.depth./cZ_data);
% depth_data = interp1(time_data, core_synthetic.depth, 0.5*radar.TWTT);

% % Calculate the residuals and variance in depth estimates between the 
% % measured density in the synthetic core and the modeled density
% depth_res = depth_data - radar.depth;
% radar.depth_var = movvar(depth_res, round(length(depth_res)/5));

end