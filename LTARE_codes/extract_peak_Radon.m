function [GRad_peaksr,GRad_peaksc]=extract_peak_Radon(GRad,psm,maxpeaks)
%calling: [GRad_peaksr,GRad_peaksc]=extract_peak_Radon(GRad,psm)
%input parameters: GRad  : the radon transform domain
%                  psm   : the percentage spot maximum, typically in [.1, .3] 
%                  
%output parameters: GRad_peaksr,GRad_peaksc: contain all spot points >= psm
%of max relatively to rows or colomns, Here we use only relatively to
%columns (GRad_peaksc)
%                   
%Description: extract_peak_Radon find all spot points >= psm of max


[nr,nc] = size(GRad);
GRad_peaksc = zeros(nr,nc);
GRad_peaksr = zeros(nr,nc);

vmax = max(max(GRad));
[r,c]=find(GRad == vmax);

Tmax = psm*vmax;


    if(size(GRad(:,c),2) >= 1)       
            %[locsc] = peakfinder(GRad(:,c(1)), Tmax); %ok 
            [locsc] = peakfinder(GRad(:,c(1)), (max(GRad(:,c(1)))-min(GRad(:,c(1))))/maxpeaks, psm*max(GRad(:,c(1))), 1);
            %[locsc] = peakfinder(GRad(:,c), (max(GRad(:,c))-min(GRad(:,c)))/20, psm*max(GRad(:,c)), 1);
            
            %[locsr] = peakfinder(GRad(r,:),Tmax);%, 1);
            %[locs2] = peakfinder(GRad(:,c), Tmax2);%, 1);
            %[rlocs,clocs] = find(GRad(:,:) >= 0.95*vmax),
           
            if(~isempty(locsc))
                for k=1:size(c,1)
                 GRad_peaksc(locsc,c(k)) = 100; %flag peaks located at all locs to be trace
                %GRad_peaks(r,locs) = 100;
                end
            end

    end
    
    if(size(GRad(r,:),1) >= 1)       
            %[locsc] = peakfinder(GRad(:,c), Tmax); %ok 
            %[locsr] = peakfinder(GRad(r(1),:),Tmax);%, 1);
            [locsr] = peakfinder(GRad(:,r(1)), (max(GRad(:,r(1)))-min(GRad(:,r(1))))/maxpeaks, psm*max(GRad(:,r(1))), 1);
            
            %[locs2] = peakfinder(GRad(:,c), Tmax2);%, 1);
            %[rlocs,clocs] = find(GRad(:,:) >= 0.95*vmax),
            
            if(~isempty(locsr))
                %GRad_peaks (locsc,c) = 100; %flag peaks located at all locs to be trace
                GRad_peaksr(r(1),locsr) = 100;

            end           
        
    end

 
end