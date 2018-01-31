function [ima_layers, over_amp, stats, peaks, ppick, pindex,header_ppick final_pick_ave, ppick_ave,pindex_ave, header_ppick_ave]=...
    LTARE_ppick(data, header, bsize,maxlayer,L0,psm,maxpeaks,delta,cols_to_ave)
%%
bsize=50; maxlayer=50; L0=3; psm=0.21; maxpeaks=16.5; delta=13; cols_to_ave=100; dist=5; min_height=0.025;

%The Layer Tracking Algorithm from Radar Echograms: LTARE
%**With the addition of the Peak Pick script to pick out the most accurate layers and their index values over each assigned column
%
%Input parameters:
%(for LTARE)        file_in,file_out: input file and output file respectivelyd
%
%                   bsize: block size, typically in [50, 100] pixels, set default to 50
%
%                   maxlayer: max layer segments to extract => [10,50], set default to 50 if more than 10 layers and 10 if less than 5
%
%                   L0: intergration length, typically in [1, 5] pixels, set default to 5
%
%                   psm: percentage of spot maximum to consider, typically in [0.1, 0.3], set default to 0.1
%
%                   maxpeaks: maximum nber of peaks in the radon space => [4,100], set default to 20
%
%                   delta: minimum distance between points to be considered as belonging to the same layer ==> [5,30]: default set 5
%
%(for Peak Pick)    col_to_ave: the width (in blocks) of columns to calculate the average for independantly
%
%
%
%Output parameters:
%(from LTARE)       ima_layers:  picked layers after radon transform,
%           nearest neight tenertopation and cleaning short segments
%
%                   over_amp  : echogram overlaped by the result layers
%                   from ima_layers
%
%
%
%(from Peak Pick)
%
%
%Description: the LTARE computes the radon transform (forward and inverse) by
%blocks of bsize X bsize, using L0 integration length, and psm = percentage
%of spot maximum to consider.
%The LTARE has been tested with both ku-band radar and radar depth sounder
%echograms
%==========================================================================

% assign data and convert to power

ima_in = 20*log10(data);
%dist=dist_KBrunt(header(3,:), header(4,:)); %set distance for images
%xt=1:floor(length(header(3,:))/5):length(header(3,:));
%xt=[xt, length(header(3,:))];


%=========================================================================


%+++++++++++++Define variables++++++++++++++++++++++++++++++++++++++++++++
[nr,nc] = size(ima_in);
pnr = ceil(nr/bsize)*bsize;
pnc = ceil(nc/bsize)*bsize;
ynr = floor((pnr+1)/2);
xnc = floor((pnc+1)/2);

raw_layers = zeros(pnr,pnc); %raw layers
ima_layers =zeros(pnr,pnc); % final results after smoothing raw layers
ima_layers2 =zeros(pnr,pnc); % final results after smoothing using spline
ima_layer = zeros(pnr,pnc); % individual layer extraction
block = zeros(bsize,bsize); %block of bsizeXbsize window for processing the RT
over_amp = zeros(pnr,pnc);
over_amp(1:nr,1:nc) = ima_in(1:nr,1:nc);%overlaping original image with resulting layer extraction

%delta between y values must be
%<= 5 to be considered as the same layer ku-band
%<= 2 for radar depth sounder
%delta = 2; %for radar depth sonder
%delta = 5; %5 for ku-band radar
%=========================================================================

%++++++++++++++++++Angular directions for ku-band radar++++++++++++++++++++
%theta6 = -pi/6:pi/(3*bsize):pi/6-pi/(3*bsize); %ok
theta5 = -pi/5:pi/(2.5*bsize):pi/5-pi/(2.5*bsize); %ok+ for internal layers
%theta1 = -pi/2:pi/(1*bsize):pi/2-pi/(1*bsize); %ok

%++++++++++++++++++Angular directions for radar depth sounder++++++++++++++
%theta64 = -pi/64:pi/(32*bsize):pi/64-pi/(32*bsize); %ok for surface/bottom

%Set the correct angular direction according to the data type
theta = theta5;
%theta = theta64;
%==========================================================================

blck = 0;%count the number of blocks
for r=1:bsize:pnr
    for c=1:bsize:pnc
        %get a block of data
        block(:,:) = over_amp(r:r+bsize-1,c:c+bsize-1);
        
        %forward radon transform
        [GRad]=local_GRad(block,L0,theta);
        
        %backward radon transform
        [Orig_block]=local_iGRad(GRad,block,L0,psm,theta,maxpeaks);
        
        %place the result in layer map
        raw_layers(r:r+bsize-1,c:c+bsize-1) = Orig_block;
        
        blck = blck+1, %prints to screen for progress
    end
end


%++++++++++++++++++++++++++ smoothing ++++++++++++++++++++++++++++++++++++

%%%collect y points
lpy = zeros(pnr,pnc); %contains all y values
lpx = zeros(pnr,pnc); %contains all x values
for l=1:maxlayer % max number of layers set
    px = []; py = [];
    for x=1:pnc
        for y=1:pnr
            if(raw_layers(y,x) ~= 0 && isempty(py)) %there is a point (get the first y)
                py = [py y];
                px = [px x];
                raw_layers(y,x) = 0; %clear the current point
                break;
            end
            
            if(raw_layers(y,x) ~= 0 && ~isempty(py)) % get a next y
                if(abs(py(size(py,2)) - y) <= delta) %delta between y must be <= 2 for Deep surface
                    %if(abs(py(size(py,2)) - y) <= 5) %delta between y must be
                    %<= 5 to consider the same layer ku-band
                    py = [py y];
                    px = [px x];
                    raw_layers(y,x) = 0; %clear the current point
                    break;
                else
                    
                    continue;
                end
            end
        end%for y
    end%for x
    
    lpy(l,1:size(py,2)) = py(1:size(py,2)); %get all y for a given layer
    lpx(l,1:size(px,2)) = px(1:size(px,2)); %get all x for a given layer
end%for l

%% cleaning to keep correct points
for l=1:maxlayer
    l, %prints to screen for progress
    cpx = lpx(l,:);
    cpy = lpy(l,:);
    
    [ppy,ppx] = clean_zeros_py(cpy,cpx);
    
    %eliminate short segments less the 50 % of columns
    %compile the result layer tracking
    %make the overlaping image with layers
    if(~isempty(ppy) && ~isempty(ppx))
        
        ssx = find(ppx);
        if(~isempty(ssx) && length(ssx) >= 0.5*nc)
            
            
            for k=1:length(ppx)
                if(ppy(k) > 0 && ppx(k) > 0)
                    ima_layer(ppy(k),ppx(k)) = 1;
                    
                    over_amp(ppy(k),ppx(k)) = max(max(over_amp));
                end%if
            end%for
            
            ima_layers(1:pnr,1:pnc) = ima_layers(1:pnr,1:pnc) + ima_layer(1:pnr,1:pnc);
            
        end%if
        
    end%if
    
end%for l

%% ***************** layer relabelling *****************
%*****************************************************
raw_layers = ima_layers;
ima_layers = zeros(pnr,pnc);
ima_layer = zeros(pnr,pnc);

%%%collect y points
lpy = zeros(pnr,pnc); %contains all y values
lpx = zeros(pnr,pnc); %contains all x values
for l=1:maxlayer % max number of layers set, can be changed
    l,
    px = []; py = [];
    for x=1:pnc
        for y=1:pnr
            if(raw_layers(y,x) ~= 0 && isempty(py)) %there is a point (get the first y)
                py = [py y];
                px = [px x];
                raw_layers(y,x) = 0; %clear the current point
                break;
            end
            
            if(raw_layers(y,x) ~= 0 && ~isempty(py)) % get a next y
                if(abs(py(size(py,2)) - y) <= delta) %delta between y must be <= 2 for Deep surface
                    %if(abs(py(size(py,2)) - y) <= 5) %delta between y must be
                    %<= 5 to consider the same layer ku-band
                    py = [py y];
                    px = [px x];
                    
                    raw_layers(y,x) = 0; %clear the current point
                    break;
                else
                    
                    continue;
                end
            end
        end%for y
    end%for x
    
    lpy(l,1:size(py,2)) = py(1:size(py,2)); %get all y for a given layer
    lpx(l,1:size(px,2)) = px(1:size(px,2)); %get all x for a given layer
end%for l

%% 2nd cleaning to keep correct points
for l=maxlayer:-1:1
    l,
    cpx = lpx(l,:);
    cpy = lpy(l,:);
    
    ppx = cpx;
    ppy = cpy;
    
    if(~isempty(ppy) && ~isempty(ppx))
        ssx = find(ppx);
        if(~isempty(ssx) && length(ssx) > 0.5*nc)
            
            
            for k=1:length(ppx)
                if(ppy(k) > 0 && ppx(k) > 0)
                    ima_layer(ppy(k),ppx(k)) = 1;
                    
                    over_amp(ppy(k),ppx(k)) = max(max(over_amp));
                end%if
            end%for
            
            ima_layers(1:pnr,1:pnc) = ima_layers(1:pnr,1:pnc) + ima_layer(1:pnr,1:pnc);
            
        end%if
        %toc,
        %pause
    end%if
    
    
    
end%for l



%************************************


%display the overlaping image with layers
figure
imagesc(over_amp)

%dispay the result layer tracking
figure
imagesc(ima_layers)


%end of LTARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Peak Picker

% specifies the beginning and end to each fraction of the data being evaluated seperately
start=1:cols_to_ave:size(ima_in,2);
stop=cols_to_ave:cols_to_ave:size(ima_in,2);

% creates cells for picks and associated index values
ppick=cell(length(start),1);
pindex=cell(length(start),1);

% determines midpt for each column
midpt=[cols_to_ave/2:cols_to_ave:size(ima_in,2)];

% returns statistics on the picked layers as well as the peaks where the
% rows were picked more often on both layers and layers2.  Layers are from
% the original LTARE and layers2 has been nearest neighbor interpolated

for i=1:maxlayer
    if ~isempty(find(ima_layers==i))
        [row,col]=find(ima_layers==i);
        stats(:,i)=[mean(row), max(row), min(row), std(row), (length(row)./size(ima_layers,2))];
    elseif isempty(find(ima_layers==i))
    end
end


for i=1:length(start)
    peaks(i,:)=mean(ima_layers(:,start(i):stop(i)),2);
    ind=find(peaks(i,:)>0);
    for j=1:(length(ind)-2)
        if peaks(i,ind(j))==peaks(i,ind(j)+1) | peaks(i,ind(j))==peaks(i,ind(j)+2)
            peaks(i,ind(j))=peaks(i,ind(j))+0.5;
        end%we can not have  a peak with two of the same numbers so add 0.5 if the number is the same as neighbor
    end
   
    for j=1:(length(ind)-2)
        if peaks(i,ind(j))==peaks(i,ind(j)+1)
            peaks(i,ind(j))=peaks(i,ind(j))+0.5;
        end%we can not have  a peak with two of the same numbers so add 0.5 if the number is the same as neighbor
    end
    [ppick_in,pindex_in]=findpeaks(peaks(i,:),'MINPEAKDISTANCE',delta,'MINPEAKHEIGHT',0.025);
    %[ppick_in,pindex_in]=peakfinder(peaks(i,:));
    ppick{i}=ppick_in;
    pindex{i}=pindex_in;
end

%get the header information for the points of ppick The header is of the
%exact point
header_ppick(1,:)=header(1,midpt);
header_ppick(2,:)=header(2,midpt);
header_ppick(3,:)=header(3,midpt);
header_ppick(4,:)=header(4,midpt);
header_ppick(5,:)=header(5,midpt);
header_ppick(6,:)=header(6,midpt);
header_ppick(7,:)=header(7,midpt);


% plots the peaks and the picked peaks for final i with layers in blue and
% layers2 (smoothed) in red

figure
plot(peaks(i,:))
hold on
plot(pindex{i},ppick{i},'x')
%plot(ppick{i},pindex{i},'x')
hold off

% plots the picked peaks over the data for each column, along with the layers LTARE.m outputs

figure
imagesc(over_amp)
hold on
for i=1:length(start)
    plot(midpt(i),[pindex{i}],'kx') % plots picked peaks at midpt of each column
    vline(start(i),'k') %creates verticle lines seperating columns
end
hold off

%% Peak Picker 2 for the entire distance (largest scale)

pick_all=pindex{1}; %load in all of the picks
for i=2:length(pindex)
    pick_all=[pick_all, pindex{i}];
end

%hist(pick_all,size(ima_layers,2))
edges=[1:size(ima_layers,1)]; %define the bins of the histogram
[final_pick_ave,bin] = histc(pick_all,edges); %get the count per bin

%pick peaks
ind=find(final_pick_ave>0);
for j=1:(length(ind)-2)
    if final_pick_ave(ind(j))==final_pick_ave(ind(j)+1)|final_pick_ave(ind(j))==final_pick_ave(ind(j)+2)|final_pick_ave(ind(j))==final_pick_ave(ind(j)+3)
        final_pick_ave(ind(j))=final_pick_ave(ind(j))+0.5;
    end
end

for j=1:(length(ind)-2)
    if final_pick_ave(ind(j))==final_pick_ave(ind(j)+1)|final_pick_ave(ind(j))==final_pick_ave(ind(j)+2)
        final_pick_ave(ind(j))=final_pick_ave(ind(j))+0.5;
    end
end

for j=1:(length(ind)-1)
    if final_pick_ave(ind(j))==final_pick_ave(ind(j)+1)
        final_pick_ave(ind(j))=final_pick_ave(ind(j))+0.5;
    end
end

[ppick_ave_in,pindex_ave_in]=findpeaks(final_pick_ave,'MINPEAKDISTANCE',delta,'MINPEAKHEIGHT',1.99);
ppick_ave=ppick_ave_in;
pindex_ave=pindex_ave_in;

%plot averages and peak picks
figure
plot(final_pick_ave)
hold on
plot(pindex_ave,ppick_ave,'x')
hold off
xlabel('Bin Number')
ylabel('Number of Picks')

%plot horizontal lines over data
figure
imagesc(over_amp)
hold on
hline(pindex_ave,'b')
hold off
set(gca,'XTick',xt,'XTickLabel',round2(dist(xt),.01))
xlabel('Distance (km)')
ylabel('Bin Number')

%get header information for the final pindex_ave data.  This is the central
%point not an average.
[nr,nc]=size(header_ppick);

header_ppick_ave=header_ppick(:,floor(nc/2));

end