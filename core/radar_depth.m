% This function provides an estimated depth scale for radar images based
% off of depth-density relations from nearby firn cores

function [radar] = radar_depth(radar, cores)

% Firn core data must be provided in the proper format (i.e. as outputted
% by "import_cores.m"

% Indices for composite core locations (every 10 km along echogram)
comp_idx = 1:round(10000/mean(diff(radar.dist))):size(radar.data_out, 2);
% rho_coeff = zeros(4, length(comp_idx));
rho_coeff = zeros(3, length(comp_idx));
rho_var = zeros(3, length(comp_idx));
depths = zeros(size(radar.data_out, 1), length(comp_idx));

for i = 1:length(comp_idx)
    
    % Generate a core composite based on the average of spatially-weighted
    % nearby cores
    [core_composite] = rho_spW(radar.Easting(comp_idx(i)), ...
        radar.Northing(comp_idx(i)), cores);
    rho_var(:,i) = core_composite.rho_var;

    
    
    rho_ice = 0.917;
    cutoff = 401;
    X = core_composite.depth;
    Y = log(core_composite.rho./(rho_ice-core_composite.rho));
    
    p0 = polyfit(X(1:cutoff+50), Y(1:cutoff+50), 1);
    p1 = polyfit(X(cutoff-50:end), Y(cutoff-50:end), 1);
    x_idx = round((p0(2)-p1(2))/(p1(1)-p0(1))/0.02);
    
    depth0 = 0:0.02:core_composite.depth(x_idx);
    depth1 = core_composite.depth(x_idx+1):0.02:150;
    rho0 = (exp(p0(1)*depth0+p0(2))./(1+exp(p0(1)*depth0+p0(2))))*rho_ice;
    rho1 = (exp(p1(1)*depth1+p1(2))./(1+exp(p1(1)*depth1+p1(2))))*rho_ice;
    
    depth_mod = [depth0 depth1]';
    rho_mod = [rho0 rho1]';
    
    
%     g = @(a, b, c, x) a*(x).^b + c;
%     rho_fit = fit(core_composite.depth(25:end), core_composite.rho(25:end), g,...
%         'StartPoint', [1 0.5 mean(core_composite.rho(1:25))]);
%     depth_mod = (0:0.02:150)';
%     rho_coeff(:,i) = coeffvalues(rho_fit);
%     rho_mod = rho_coeff(1,i)*depth_mod.^rho_coeff(2,i) + rho_coeff(3,i);
    
    % Diagnostic figures for density modeling
    figure
    hold on
    plot(core_composite.rho, core_composite.depth)
    plot(rho_mod, depth_mod)
    ylim([0 core_composite.depth(end)])
    ylabel('Depth (m)')
    xlabel('Density (g/cm^3)')
    set(gca, 'Ydir', 'reverse')
    legend('Composite firn core', 'Density model')
    hold off
    
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
    radar.TWTT = (0:2*radar.time_trace(1):2*radar.time_trace(1)*...
        (size(radar.data_out, 1)-1))';
    
    % Convert recorded TWT time to depth by interpolating measured values to
    % the modeled time-depth relationship
    depths(:,i) = interp1(time_mod, depth_mod, 0.5*radar.TWTT);
end

% Interpolate between data points for core composites with more than 1
% component
if length(comp_idx) > 1
    radar.rho_coeff = (interp1(comp_idx, rho_coeff', ...
        1:size(radar.data_out, 2), 'linear', 'extrap'))';
    radar.rho_var = (interp1(comp_idx, rho_var', ...
        1:size(radar.data_out, 2), 'linear', 'extrap'))';
    radar.depth = (interp1(comp_idx, depths', ...
        1:size(radar.data_out, 2), 'linear', 'extrap'))';
else
    radar.rho_coeff = rho_coeff;
    radar.rho_var = rho_var;
    radar.depth = depths;
end

end
