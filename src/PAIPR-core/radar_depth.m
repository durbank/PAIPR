% This function provides an estimated depth scale for radar images based
% off of depth-density relations from nearby firn cores

function [radar] = radar_depth(radar, rho_data)

% Firn core data must be provided in the proper format (i.e. as outputted
% by "import_cores.m"

% Indices for composite core locations (every 10 km along echogram)
depths = zeros(size(radar.data_stack));

for i = 1:size(radar.data_stack,2)
    
    % Extract depth and predicted mean density from saved rho data
    depth_mod = rho_data.Data{i}.Depth;
    rho_mod = rho_data.Data{i}.pred_mean;

    % Calculate the real part of the relative permittivity of firn
    e_mod = (1 + 0.845*rho_mod).^2;  % Kovacs
    %     e_mod2 = ((rho_mod/rho_ice)*...
    %         (emiss_ice^(1/3)-1) + 1).^3; % Looyenga, 1965
    
    % Calculate radar propagation speed with depth from emissivity
    c0 = 2.9979E8;  % speed of light in a vacuum (m/s)
    cZ_mod = c0./sqrt(e_mod);   % Speed of light with depth in firn (m/s)
    
    % Calculate modeled one way travel time (seconds)
%     time_discrete = mean(diff(depth_mod))./cZ_mod;
%     time_mod = cumsum([0; time_discrete(1:end-1)]);
    time_mod = cumsum([0; diff(depth_mod)]./cZ_mod);
    
    % Calculate cumulative two-way travel time for radar data
    TWTT = (0:2*radar.time_trace(1):2*radar.time_trace(1)*...
        (size(radar.data_stack, 1)-1))';
    
    % Convert recorded TWT time to depth by interpolating measured values 
    % to the modeled time-depth relationship
    depths(:,i) = interp1(time_mod, depth_mod, 0.5*TWTT, ...
        'pchip', 'extrap');
end

% Add depth variable to radar structure
radar.depth = depths;

end
