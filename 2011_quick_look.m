% to run this code declare output directory on line 16 and on line 12 input the
% filenames that you would like to take a quick look at in the field.

% Ben Panzer
% William Blake
% Carl Leuschen
%adapted by Lora Koenig for WAIS traverse 2011.

clear all; close all; clc;
format compact; format short;

filename=struct('name', {'kuband_20110923_15325088_0002.dat',...
    'kuband_20110923_15314382_0001.dat'});

% Output Directory
out_dir = '/icebridgedata/lorak/test/';

% Filebase for the data

pulse_length = 250e-6;
kubandwidth = 5.0e9; %ku_band
snowbandwidth=6.0e9; %snow radar

% Title for echogram
info_str   = 'Quick look radar WAIS 2011';

%XXXXXXXXXXXXXXXXXXXX NOT NECESSARY TO MODIFY FOR ARCTIC DATA XXXXXXXXXXXXXXXXXXXXXXX%

% Sampling frequency of the data acquisition system
fs = 62.5e6;

% Length of the FFT %lk or record length  this has to be the length of the
% records recorded by the radar or it will give the wrong pixel size.
%FFT_len = floor((pulse_length*fs)/1000)*1000;
FFT_len=16000;


% Real part of the dielectric constant of dry snow (currently set for air)

rho=0.320;
eps_snow = 1+1.5995*rho+1.861*rho^3; %from Matzler
% filenames that you would like to take a quick look at in the field.

% Ben Panzer
% William Blake
% Carl Leuschen
%adapted by Lora Koenig for WAIS traverse 2011.

clear all; close all; clc;
format compact; format short;

filename=struct('name', {'kuband_20110923_15325088_0002.dat',...
    'kuband_20110923_15314382_0001.dat'});

% Output Directory
out_dir = '/icebridgedata/lorak/test/';

% Filebase for the data

pulse_length = 250e-6;
kubandwidth = 5.0e9; %ku_band
snowbandwidth=6.0e9; %snow radar

% Title for echogram
info_str   = 'Quick look radar WAIS 2011';

%XXXXXXXXXXXXXXXXXXXX NOT NECESSARY TO MODIFY FOR ARCTIC DATA XXXXXXXXXXXXXXXXXXXXXXX%

% Sampling frequency of the data acquisition system
fs = 62.5e6;

% Length of the FFT %lk or record length  this has to be the length of the
% records recorded by the radar or it will give the wrong pixel size.
%FFT_len = floor((pulse_length*fs)/1000)*1000;
FFT_len=16000;


% Real part of the dielectric constant of dry snow (currently set for air)

rho=0.320;
eps_snow = 1+1.5995*rho+1.861*rho^3; %from Matzler

% Number of presums is based on the unfocused synthetic aperture length
presums    = 4;

%XXXXXXXXXXXXXXXXXXXXXX END OF PARAMETERS XXXXXXXXXXXXXXXXXXXXXXXXXXX%

% Pixel size is a function of the length of FFT, pulse length, bandwidth
% and the assumed dielectric constant of dry snow
c   = 299792458;
ku_pixel_size_snow = (fs*pulse_length*c)/(2*FFT_len*kubandwidth*sqrt(eps_snow));
ku_pixel_size_air = (fs*pulse_length*c)/(2*FFT_len*kubandwidth);
ku_delta_t = (fs*pulse_length)/(2*FFT_len*kubandwidth); %use this as pixel size generic

snow_pixel_size_snow = (fs*pulse_length*c)/(2*FFT_len*snowbandwidth*sqrt(eps_snow));
snow_pixel_size_air = (fs*pulse_length*c)/(2*FFT_len*snowbandwidth);
snow_delta_t = (fs*pulse_length)/(2*FFT_len*snowbandwidth); %use this as pixel size generic

for i1 = 1:length(filename)
    tic
    fid = fopen(filename(i1).name,'r','ieee-be');
    deadbeef = hex2dec('deadbeef');
    
    while 1
        tmp_search     = fread(fid,80000,'uint8');
        tmp_search_len = length(tmp_search);
        tmp_search     = tmp_search(1:(floor(tmp_search_len/4)*4));
        tmp            = tmp_search(1:tmp_search_len-3)*2^24+tmp_search(2:tmp_search_len-2)*2^16+ ...
            tmp_search(3:tmp_search_len-1)*2^8+tmp_search(4:tmp_search_len);
        if isempty(tmp_search)
            error('No Records found in %s',just_fn);
        end
        deadbeef_idx = find(tmp == deadbeef);
        if length(deadbeef_idx) >=2
            rec_len = (deadbeef_idx(2) - deadbeef_idx(1))/2;
        end
        if ~isempty(deadbeef_idx)
            fseek(fid,-tmp_search_len,'cof');
            fseek(fid,(deadbeef_idx(1)-1)*1,'cof');
            break;
        end
        fseek(fid,-4,'cof');
    end
    index = ftell(fid);
    clear tmp tmp_search tmp_search_len
    
    % Windowing for fast and slow time
    fast_time_win = hanning(FFT_len/2); %LK fast time is the sampling rate
    %divide by two for the two radars
    slow_time_win = hanning(presums); %LK slow time is the PRF May want
    %to remove
    
    % Resampling the data as samples are flipped
    %%LK are the data I have  flipped? Yes Per Ben and Aqsa This is a firm ware
    %%issue with the radars and is fixed here.  
    rsmp = zeros(1,FFT_len);
    for i0 = 1:FFT_len/2,
        rsmp(2*i0-1) = 2*i0 + 16; %+16 is for hte header data
        rsmp(2*i0) = 2*i0-1 + 16;
    end;
    
    % Read in data
    fseek(fid,index,'bof');
    fprintf('Reading In data\n');
    rec_len=16016; %lk hard coded for snow and Ku data in the same file
    Data = fread(fid,[rec_len,inf],'uint16=>float32');
    fclose(fid);
    
    fprintf('Processing data\n');
    hdr = Data(1:16,:);
    Data = Data(rsmp,:);
    
    % Random records contain saturated triplets spaced by 3 and 2 or
    % quadruplets spaced by 2, 1, and 2.
    % Find these triplets and replace them with the mean of adjacent
    % samples
    for idx = 1:size(Data,2)
        bad_idx = find(Data(:,idx) < 4000);
        if (length(bad_idx) == 3) && all(diff(bad_idx) == [3; 2])
            if bad_idx(1) > 1
                Data(bad_idx(1),idx) = (Data(bad_idx(1)-1,idx)+Data(bad_idx(1)+1,idx))/2;
            else
                Data(bad_idx(1),idx) = Data(bad_idx(1)+1,idx);
            end
            Data(bad_idx(2),idx) = (Data(bad_idx(2)-2,idx)+Data(bad_idx(2)+1,idx))/2;
            Data(bad_idx(2)-1,idx) = (Data(bad_idx(2)-2,idx)+Data(bad_idx(2)+1,idx))/2;
            if bad_idx(3) < size(Data,1)
                Data(bad_idx(3),idx) = (Data(bad_idx(3)-1,idx)+Data(bad_idx(3)+1,idx))/2;
            else
                Data(bad_idx(3),idx) = Data(bad_idx(3)-1,idx);
            end
        elseif (length(bad_idx) == 4) && all(diff(bad_idx) == [2; 1; 2])
            if bad_idx(1) > 1
                Data(bad_idx(1),idx) = (Data(bad_idx(1)-1,idx)+Data(bad_idx(1)+1,idx))/2;
            else
                Data(bad_idx(1),idx) = Data(bad_idx(1)+1,idx);
            end
            Data(bad_idx(3),idx) = (Data(bad_idx(3)-2,idx)+Data(bad_idx(3)+1,idx))/2;
            Data(bad_idx(2),idx) = (Data(bad_idx(2)-1,idx)+Data(bad_idx(2)+2,idx))/2;
            if bad_idx(4) < size(Data,1)
                Data(bad_idx(4),idx) = (Data(bad_idx(4)-1,idx)+Data(bad_idx(4)+1,idx))/2;
            else
                Data(bad_idx(4),idx) = Data(bad_idx(4)-1,idx);
            end
        end
    end
    
    %seperate ku and snow radar data at this point and process
    
    kuData=Data(1:8000,:);
    snowData=Data(8001:16000,:);
    clear Data
    
    %from here are everything is processed in duplicate
    
    % Subtract the mean from the data
    for st_idx = 1:size(kuData,2)
        kuData(:,st_idx) = kuData(:,st_idx)-mean(kuData(:,st_idx));
    end
    
     % Subtract the mean from the data
    for st_idx = 1:size(snowData,2)
        snowData(:,st_idx) = snowData(:,st_idx)-mean(snowData(:,st_idx));
    end
    
    % Perform Coherent Integrations ku
    siz = size(kuData);
    new_len      = floor(siz(2)/presums);
    nearest_len  = floor(siz(2)/presums)*presums;
    kuData         = reshape(kuData(:,1:nearest_len), ...
        [siz(1) presums new_len]);
    tmp_kuData      = zeros(size(kuData,1),size(kuData,3),'single');
    slow_time_win = repmat(slow_time_win.',[size(kuData,1),1,1]);
    for st_idx = 1:size(kuData,3);
        kuData(:,:,st_idx)   = kuData(:,:,st_idx).*slow_time_win;
        tmp_kuData(:,st_idx) = mean(kuData(:,:,st_idx),2);
    end
    kuData = double(tmp_kuData);
    clear tmp_kuData slow_time_win
    
    slow_time_win = hanning(presums); %LK slow time is the PRF May want
    %to remove
    
    % Perform Coherent Integrations snow
    siz = size(snowData);
    new_len      = floor(siz(2)/presums);
    nearest_len  = floor(siz(2)/presums)*presums;
    snowData         = reshape(snowData(:,1:nearest_len), ...
        [siz(1) presums new_len]);
    tmp_snowData      = zeros(size(snowData,1),size(snowData,3),'single');
    slow_time_win = repmat(slow_time_win.',[size(snowData,1),1,1]);
    for st_idx = 1:size(snowData,3);
        snowData(:,:,st_idx)   = snowData(:,:,st_idx).*slow_time_win;
        tmp_snowData(:,st_idx) = mean(snowData(:,:,st_idx),2);
    end
    snowData = double(tmp_snowData);
    clear tmp_snowData slow_time_win 
    
     % Window data in fast time ku
    for st_idx = 1:size(kuData,2)
        kuData(:,st_idx) = kuData(:,st_idx).*fast_time_win;
    end
    
     % Window data in fast time snow
    for st_idx = 1:size(snowData,2)
        snowData(:,st_idx) = snowData(:,st_idx).*fast_time_win;
    end
   
    
   %fft ku data 
    kuData     = fft(kuData);
    fft_idxs = 1:(floor(size(kuData,1)/2));
    kuData     = kuData(fft_idxs,:);
    kuData     = abs(kuData);
    
   %fft snow data
    snowData     = fft(snowData);
    fft_idxs = 1:(floor(size(snowData,1)/2));
    snowData     = snowData(fft_idxs,:);
    snowData     = abs(snowData);
    
%       
%         % Tracking the surface return and eliminating outliers to set the data
%         % start and stop range
%    
% 
%     %LK for the ground based radars we know where the surface is and the the
%     %penetration depth so we can set the istart and istop values.  The
%     istart and istop values are set for ku and snow seperately to
%     seperater them in the file.  The file records the Kuband first and
%     then the snow.
   

    istart = 1; %allow us to see the antenna height at 170 or 171 is the snow surface     
    istop = 4000; %this is 40 meters into the snow assuming rho=320. and is below the penetration depth.
    
    kuData = kuData(istart:istop,:);
    snowData            = snowData(istart:istop,:);
    kuDepth           = ((0:istop-istart))*ku_pixel_size_snow;
    snowDepth = ((0:istop-istart))*snow_pixel_size_snow;
    

    fprintf('Creating images\n');
    
    % Plot Normal
    hf = figure;
    imagesc([],kuDepth,20.*log10(kuData));
    caxis([40 150]);
    hold on;
    colormap(1-bone);
    
 
    img_title = sprintf('Quick Look Kuband %s',filename(i1).name);
    title(img_title,'FontWeight','Bold');
    ylabel('Depth [m]');
    xlabel('Distance [column number]');
    picfilename = sprintf('%sKuband_FFT_image.%02d.%04d',out_dir,filename(i1).name,i1);
    print('-djpeg','-r300',[picfilename '.jpg']);
    close(hf);

    
        % Plot Normal
    hf = figure;
    imagesc([],snowDepth,20.*log10((snowData)));
    caxis([40 150]);
    hold on;
    colormap(1-bone);
    
 
    img_title = sprintf('Quick Look Snow %s',filename(i1).name);
    title(img_title,'FontWeight','Bold');
    ylabel('Depth [m]');
    xlabel('Distance [column number]');
    picfilename = sprintf('%sSnow_FFT_image.%02d.%04d',out_dir,filename(i1).name,i1);
    print('-djpeg','-r300',[picfilename '.jpg']);
    close(hf);
    
    toc
end
% Number of presums is based on the unfocused synthetic aperture length
presums    = 4;

%XXXXXXXXXXXXXXXXXXXXXX END OF PARAMETERS XXXXXXXXXXXXXXXXXXXXXXXXXXX%

% Pixel size is a function of the length of FFT, pulse length, bandwidth
% and the assumed dielectric constant of dry snow
c   = 299792458;
ku_pixel_size_snow = (fs*pulse_length*c)/(2*FFT_len*kubandwidth*sqrt(eps_snow));
ku_pixel_size_air = (fs*pulse_length*c)/(2*FFT_len*kubandwidth);
ku_delta_t = (fs*pulse_length)/(2*FFT_len*kubandwidth); %use this as pixel size generic

snow_pixel_size_snow = (fs*pulse_length*c)/(2*FFT_len*snowbandwidth*sqrt(eps_snow));
snow_pixel_size_air = (fs*pulse_length*c)/(2*FFT_len*snowbandwidth);
snow_delta_t = (fs*pulse_length)/(2*FFT_len*snowbandwidth); %use this as pixel size generic

for i1 = 1:length(filename)
    tic
    fid = fopen(filename(i1).name,'r','ieee-be');
    deadbeef = hex2dec('deadbeef');
    
    while 1
        tmp_search     = fread(fid,80000,'uint8');
        tmp_search_len = length(tmp_search);
        tmp_search     = tmp_search(1:(floor(tmp_search_len/4)*4));
        tmp            = tmp_search(1:tmp_search_len-3)*2^24+tmp_search(2:tmp_search_len-2)*2^16+ ...
            tmp_search(3:tmp_search_len-1)*2^8+tmp_search(4:tmp_search_len);
        if isempty(tmp_search)
            error('No Records found in %s',just_fn);
        end
        deadbeef_idx = find(tmp == deadbeef);
        if length(deadbeef_idx) >=2
            rec_len = (deadbeef_idx(2) - deadbeef_idx(1))/2;
        end
        if ~isempty(deadbeef_idx)
            fseek(fid,-tmp_search_len,'cof');
            fseek(fid,(deadbeef_idx(1)-1)*1,'cof');
            break;
        end
        fseek(fid,-4,'cof');
    end
    index = ftell(fid);
    clear tmp tmp_search tmp_search_len
    
    % Windowing for fast and slow time
    fast_time_win = hanning(FFT_len/2); %LK fast time is the sampling rate
    %divide by two for the two radars
    slow_time_win = hanning(presums); %LK slow time is the PRF May want
    %to remove
    
    % Resampling the data as samples are flipped
    %%LK are the data I have  flipped? Yes Per Ben and Aqsa This is a firm ware
    %%issue with the radars and is fixed here.  
    rsmp = zeros(1,FFT_len);
    for i0 = 1:FFT_len/2,
        rsmp(2*i0-1) = 2*i0 + 16; %+16 is for hte header data
        rsmp(2*i0) = 2*i0-1 + 16;
    end;
    
    % Read in data
    fseek(fid,index,'bof');
    fprintf('Reading In data\n');
    rec_len=16016; %lk hard coded for snow and Ku data in the same file
    Data = fread(fid,[rec_len,inf],'uint16=>float32');
    fclose(fid);
    
    fprintf('Processing data\n');
    hdr = Data(1:16,:);
    Data = Data(rsmp,:);
    
    % Random records contain saturated triplets spaced by 3 and 2 or
    % quadruplets spaced by 2, 1, and 2.
    % Find these triplets and replace them with the mean of adjacent
    % samples
    for idx = 1:size(Data,2)
        bad_idx = find(Data(:,idx) < 4000);
        if (length(bad_idx) == 3) && all(diff(bad_idx) == [3; 2])
            if bad_idx(1) > 1
                Data(bad_idx(1),idx) = (Data(bad_idx(1)-1,idx)+Data(bad_idx(1)+1,idx))/2;
            else
                Data(bad_idx(1),idx) = Data(bad_idx(1)+1,idx);
            end
            Data(bad_idx(2),idx) = (Data(bad_idx(2)-2,idx)+Data(bad_idx(2)+1,idx))/2;
            Data(bad_idx(2)-1,idx) = (Data(bad_idx(2)-2,idx)+Data(bad_idx(2)+1,idx))/2;
            if bad_idx(3) < size(Data,1)
                Data(bad_idx(3),idx) = (Data(bad_idx(3)-1,idx)+Data(bad_idx(3)+1,idx))/2;
            else
                Data(bad_idx(3),idx) = Data(bad_idx(3)-1,idx);
            end
        elseif (length(bad_idx) == 4) && all(diff(bad_idx) == [2; 1; 2])
            if bad_idx(1) > 1
                Data(bad_idx(1),idx) = (Data(bad_idx(1)-1,idx)+Data(bad_idx(1)+1,idx))/2;
            else
                Data(bad_idx(1),idx) = Data(bad_idx(1)+1,idx);
            end
            Data(bad_idx(3),idx) = (Data(bad_idx(3)-2,idx)+Data(bad_idx(3)+1,idx))/2;
            Data(bad_idx(2),idx) = (Data(bad_idx(2)-1,idx)+Data(bad_idx(2)+2,idx))/2;
            if bad_idx(4) < size(Data,1)
                Data(bad_idx(4),idx) = (Data(bad_idx(4)-1,idx)+Data(bad_idx(4)+1,idx))/2;
            else
                Data(bad_idx(4),idx) = Data(bad_idx(4)-1,idx);
            end
        end
    end
    
    %seperate ku and snow radar data at this point and process
    
    kuData=Data(1:8000,:);
    snowData=Data(8001:16000,:);
    clear Data
    
    %from here are everything is processed in duplicate
    
    % Subtract the mean from the data
    for st_idx = 1:size(kuData,2)
        kuData(:,st_idx) = kuData(:,st_idx)-mean(kuData(:,st_idx));
    end
    
     % Subtract the mean from the data
    for st_idx = 1:size(snowData,2)
        snowData(:,st_idx) = snowData(:,st_idx)-mean(snowData(:,st_idx));
    end
    
    % Perform Coherent Integrations ku
    siz = size(kuData);
    new_len      = floor(siz(2)/presums);
    nearest_len  = floor(siz(2)/presums)*presums;
    kuData         = reshape(kuData(:,1:nearest_len), ...
        [siz(1) presums new_len]);
    tmp_kuData      = zeros(size(kuData,1),size(kuData,3),'single');
    slow_time_win = repmat(slow_time_win.',[size(kuData,1),1,1]);
    for st_idx = 1:size(kuData,3);
        kuData(:,:,st_idx)   = kuData(:,:,st_idx).*slow_time_win;
        tmp_kuData(:,st_idx) = mean(kuData(:,:,st_idx),2);
    end
    kuData = double(tmp_kuData);
    clear tmp_kuData slow_time_win
    
    slow_time_win = hanning(presums); %LK slow time is the PRF May want
    %to remove
    
    % Perform Coherent Integrations snow
    siz = size(snowData);
    new_len      = floor(siz(2)/presums);
    nearest_len  = floor(siz(2)/presums)*presums;
    snowData         = reshape(snowData(:,1:nearest_len), ...
        [siz(1) presums new_len]);
    tmp_snowData      = zeros(size(snowData,1),size(snowData,3),'single');
    slow_time_win = repmat(slow_time_win.',[size(snowData,1),1,1]);
    for st_idx = 1:size(snowData,3);
        snowData(:,:,st_idx)   = snowData(:,:,st_idx).*slow_time_win;
        tmp_snowData(:,st_idx) = mean(snowData(:,:,st_idx),2);
    end
    snowData = double(tmp_snowData);
    clear tmp_snowData slow_time_win 
    
     % Window data in fast time ku
    for st_idx = 1:size(kuData,2)
        kuData(:,st_idx) = kuData(:,st_idx).*fast_time_win;
    end
    
     % Window data in fast time snow
    for st_idx = 1:size(snowData,2)
        snowData(:,st_idx) = snowData(:,st_idx).*fast_time_win;
    end
   
    
   %fft ku data 
    kuData     = fft(kuData);
    fft_idxs = 1:(floor(size(kuData,1)/2));
    kuData     = kuData(fft_idxs,:);
    kuData     = abs(kuData);
    
   %fft snow data
    snowData     = fft(snowData);
    fft_idxs = 1:(floor(size(snowData,1)/2));
    snowData     = snowData(fft_idxs,:);
    snowData     = abs(snowData);
    
%       
%         % Tracking the surface return and eliminating outliers to set the data
%         % start and stop range
%    
% 
%     %LK for the ground based radars we know where the surface is and the the
%     %penetration depth so we can set the istart and istop values.  The
%     istart and istop values are set for ku and snow seperately to
%     seperater them in the file.  The file records the Kuband first and
%     then the snow.
   

    istart = 1; %allow us to see the antenna height at 170 or 171 is the snow surface     
    istop = 4000; %this is 40 meters into the snow assuming rho=320. and is below the penetration depth.
    
    kuData = kuData(istart:istop,:);
    snowData            = snowData(istart:istop,:);
    kuDepth           = ((0:istop-istart))*ku_pixel_size_snow;
    snowDepth = ((0:istop-istart))*snow_pixel_size_snow;
    

    fprintf('Creating images\n');
    
    % Plot Normal
    hf = figure;
    imagesc([],kuDepth,20.*log10(kuData));
    caxis([40 150]);
    hold on;
    colormap(1-bone);
    
 
    img_title = sprintf('Quick Look Kuband %s',filename(i1).name);
    title(img_title,'FontWeight','Bold');
    ylabel('Depth [m]');
    xlabel('Distance [column number]');
    picfilename = sprintf('%sKuband_FFT_image.%02d.%04d',out_dir,filename(i1).name,i1);
    print('-djpeg','-r300',[picfilename '.jpg']);
    close(hf);

    
        % Plot Normal
    hf = figure;
    imagesc([],snowDepth,20.*log10((snowData)));
    caxis([40 150]);
    hold on;
    colormap(1-bone);
    
 
    img_title = sprintf('Quick Look Snow %s',filename(i1).name);
    title(img_title,'FontWeight','Bold');
    ylabel('Depth [m]');
    xlabel('Distance [column number]');
    picfilename = sprintf('%sSnow_FFT_image.%02d.%04d',out_dir,filename(i1).name,i1);
    print('-djpeg','-r300',[picfilename '.jpg']);
    close(hf);
    
    toc
end

