function [CC, data_pts]=DGK_horizons(data)
%%
bsize=50; maxlayer=75; L0=3; psm=0.21; maxpeaks=16.5; delta=13; %cols_to_ave=100; dist=5; min_height=0.025;

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

% ima_in = 20*log10(data);
ima_in = data;
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




% Clipped indices of layers
hor_x = lpx(1:maxlayer,:);
hor_y = lpy(1:maxlayer,:);

% Preallocate
data_pts = zeros(size(ima_in));

% figure
% imagesc(ima_in)
% hold on
% for i = 1:maxlayer
%     plot(hor_x(i,:), hor_y(i,:), '.')
%     for j = 1:size(hor_x, 2)
%         if hor_y(i,j) <= size(ima_in, 1) && hor_y(i,j) > 0 && ...
%                 hor_x(i,j) <= size(ima_in, 2)
%             data_pts(hor_y(i,j),hor_x(i,j)) = 1;
%         end
%     end
% end
% hold off 

ima_layers = bwmorph(data_pts, 'thicken', 2);
ima_layers = bwmorph(ima_layers, 'bridge');
ima_layers = bwmorph(ima_layers, 'diag', Inf);
ima_layers = bwmorph(ima_layers, 'fill');
% ima_layers = bwmorph(ima_layers, 'thin', 2);

CC = bwconncomp(ima_layers, 4);



end
