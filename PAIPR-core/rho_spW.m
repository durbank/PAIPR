% Function to calculate the depth density profile of a radar trace from the
% spatial weighting of depth-density profiles for surrounding firn cores

function [core_composite] = rho_spW(Easting, Northing, cores)

% Range (in meters) over which SMB is thought to co-vary
range = 500000;

% Create array for radar location
loc_radar = [Easting, Northing];

% Calculate euclidean distance between cores and the current radar trace
loc_core = [(cores.Easting)' (cores.Northing)'];
loc = [loc_radar; loc_core];
mDist = squareform(pdist(loc, 'euclidean'));
vDist = mDist(1, 2:end);

% Subset all cores to only those within the defined range of
% influence/covariance
core_idx = find(vDist <= range);
cores_fd_nm = fieldnames(cores);
core_log = contains(cores_fd_nm, 'SEAT');
cores_all = cores_fd_nm(core_log);
% cores_SWM = cores_SWM(7:end);
cores_SWM = cores_all(core_idx);
vDist_SWM = vDist(core_idx);

% % Set depth of synthetic core to depth of the shallowest core used in
% % spatial weighting
% max_depth = zeros(1, numel(cores_SWM));
% for i = 1:numel(cores_SWM)
%     max_depth(i) = max(cores.(cores_SWM{i}).depth);
% end
% depth_core = (0:0.02:min(max_depth))';

% Set depth of synthetic core to depth of deepest core used in weighting
max_depth = zeros(1, numel(cores_SWM));
idx_max = zeros(1, length(cores_SWM));
for i = 1:numel(cores_SWM)
    max_depth(i) = max(cores.(cores_SWM{i}).depth);
    idx_max(i) = length(cores.(cores_SWM{i}).depth);
end
depth_core = (0:0.02:max(max_depth))';



    

% Calculate a spatial weighting matrix based on exponential inverse 
% distance (currently b = 2)
IDW = 1./vDist_SWM.^2;
SWM = IDW/sum(IDW);


IDW = zeros(length(depth_core), length(cores_SWM));
for i = 1:length(cores_SWM)  
    IDW(1:idx_max(i),i) = 1/vDist_SWM(i)^2;
end

SWM = IDW./sum(IDW,2);

% Generate a spatially weighted model of variance for the depth-density
% profile at the radar location
[rho_var_coeff] = rho_variance(cores, core_idx, SWM);

% Create data structure for the synthetic core, containing depth, density 
% age, location, and isotope values
ages = zeros(numel(depth_core), numel(cores_SWM));
rhos = zeros(numel(depth_core), numel(cores_SWM));
dDs = zeros(numel(depth_core), numel(cores_SWM));
d18Os = zeros(numel(depth_core), numel(cores_SWM));
for i = 1:numel(cores_SWM)
    age_temp = cores.(cores_SWM{i}).age;
    ages(1:length(age_temp),i) = age_temp;
%     ages(:,i) = age_temp(1:numel(depth_core));
    rho_temp = cores.(cores_SWM{i}).rho;
    rhos(1:length(rho_temp),i) = rho_temp;
%     rhos(:,i) = rho_temp(1:numel(depth_core));
    dD_temp = cores.(cores_SWM{i}).dD;
    dDs(1:length(dD_temp),i) = dD_temp;
%     dDs(:,i) = dD_temp(1:numel(depth_core));
    d18O_temp = cores.(cores_SWM{i}).d18O;
    d18Os(1:length(d18O_temp),i) = d18O_temp;
%     d18Os(:,i) = d18O_temp(1:numel(depth_core));
end
core_composite = struct('depth', depth_core, 'age', sum(SWM.*ages, 2),...
    'rho', sum(SWM.*rhos, 2), 'rho_var', rho_var_coeff, 'dD', sum(SWM.*dDs, 2), ...
    'd18O', sum(SWM.*d18Os, 2), 'Easting', loc_radar(1), 'Northing', loc_radar(2));
end