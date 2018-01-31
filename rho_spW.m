% Function to calculate the depth density profile of a radar trace from the
% spatial weighting of depth-density profiles for surrounding firn cores

function [core_synthetic] = rho_spW(radar, cores)

% Range (in meters) over which SMB is thought to covary
range = 500000;

% %%% TEMPORARY ADDITION %%%
% % This removes cores without age-depth scales from consideration (SEAT11_5
% % and SEAT11_9)
% subset = [1:9 11:13];
% cores.name = cores.name(subset);
% cores.lat = cores.lat(subset);
% cores.lon = cores.lon(subset);
% cores.elev = cores.elev(subset);
% cores.Easting = cores.Easting(subset);
% cores.Northing = cores.Northing(subset);
% %%% TEMPORARY ADDITION %%%

% Easting/Northing of the midpoint of the radar trace
loc_radar = [radar.Easting(round(numel(radar.Easting)/2)) ...
    radar.Northing(round(numel(radar.Northing)/2))];

% Calculate euclidean distance between cores and midpoint of the radar
% transect
loc_core = [(cores.Easting)' (cores.Northing)'];
loc = [loc_radar; loc_core];
mDist = squareform(pdist(loc, 'euclidean'));
vDist = mDist(1, 2:end);

% Subset all cores to only those within the defined range of
% influence/covariance
core_idx = find(vDist <= range);
cores_SWM = fieldnames(cores);
cores_SWM = cores_SWM(7:end);
% %%% TEMPORARY ADDITION %%%
% % Removes cores without age-depth scales from consideration
% cores_SWM = cores_SWM(subset);
% %%% TEMPORARY ADDITION %%%
cores_SWM = cores_SWM(core_idx);
vDist_SWM = vDist(core_idx);

% Set depth of synthetic core to depth of the shallowest core used in
% spatial weighting
max_depth = zeros(1, numel(cores_SWM));
for i = 1:numel(cores_SWM)
    max_depth(i) = max(cores.(cores_SWM{i}).depth);
end
depth_core = (0:0.02:min(max_depth))';

% % Calculate spatial weight matrix based on squared exponential kernel
% % using hyperparamters 'hp'
% hp = [1 range/exp(1)];
% cov = hp(1)^2*exp((-1/2)*vDist_SWM.^2/hp(2)^2);
% cov = cov.^2;
% SWM = cov/sum(cov);

% Calculate a spatial weighting matrix based on exponential inverse 
% distance (currently b = 2)
IDW = 1./vDist_SWM.^2;
SWM = IDW/sum(IDW);

% Create data structure for the synthetic core, containing depth, density 
% age, location, and isotope values
ages = zeros(numel(depth_core), numel(cores_SWM));
rhos = zeros(numel(depth_core), numel(cores_SWM));
dDs = zeros(numel(depth_core), numel(cores_SWM));
d18Os = zeros(numel(depth_core), numel(cores_SWM));
for i = 1:numel(cores_SWM)
    age_temp = cores.(cores_SWM{i}).age;
    ages(:,i) = age_temp(1:numel(depth_core));
    rho_temp = cores.(cores_SWM{i}).rho;
    rhos(:,i) = rho_temp(1:numel(depth_core));
    dD_temp = cores.(cores_SWM{i}).dD;
    dDs(:,i) = dD_temp(1:numel(depth_core));
    d18O_temp = cores.(cores_SWM{i}).d18O;
    d18Os(:,i) = d18O_temp(1:numel(depth_core));
end
core_synthetic = struct('depth', depth_core, 'age', sum(SWM.*ages, 2),...
    'rho', sum(SWM.*rhos, 2), 'dD', sum(SWM.*dDs, 2), ...
    'd18O', sum(SWM.*d18Os, 2), 'Easting', loc_radar(1), 'Northing', loc_radar(2));
end