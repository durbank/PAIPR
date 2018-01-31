function [Orig_ima]=local_iGRad(GRad,ima_in,L0,psm,theta,maxpeaks)
%calling : [Orig_ima]=local_iGRad(GRad,ima_in,L0,psm,theta)
%input parameters: GRad  : the radon transform domain
%                  ima_in: input image used to compute the radon transform
%                  L0    : intergration length, typically in [1, 5] pixels
%                  psm   : percentage of spot maximum to consider
%                          typicall in [0.1, 0.3]
%                  theta : angular direction range used to compute the
%                  forward radon transform 
%
%output parameters: Orig_ima: segment layers to the original space 
%                   
%
%Description: local_iGRad computes the backward radon transform  
%using GRad the forward radon transform, L0 integration length, and theta
%the angular direction range and the 
%psm = percentage of the spot max to consider



[nr,nc] = size(ima_in);
nn = min(nr,nc);
ynr = floor((nn+1)/2);
xnc = floor((nn+1)/2);

Orig_ima = zeros(nn,nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L0 = 10; %integration length
rho = -nn/2:nn/2-1;
s_rho = size(rho,2); % nn points


%theta = -pi/6:pi/(3*nn):pi/6-pi/(3*nn); %ok 
%theta = -pi/5:pi/(2.5*nn):pi/5-pi/(2.5*nn); %ok+ 
%theta = -pi/4:pi/(2*nn):pi/4-pi/(2*nn); 
%theta = -pi/3:pi/(1.5*nn):pi/3-pi/(1.5*nn); 
%theta = -pi/2:pi/nn:pi/2-pi/nn; 

s_theta = size(theta,2); %nn points
costab = cos(theta);
sintab = sin(theta);

xc = -nn/2:nn/2-1;
yr = -nn/2:nn/2-1;

%get peak locations in the radon domain
%[GRad_peaks]=extract_peak_Radon(GRad,ppc);
[GRad_peaksr,GRad_peaksc]=extract_peak_Radon(GRad,psm,maxpeaks);
%[GRad_peaks]=extract_peak_Radon(GRad,Tmax);
%[GRad_peaks]=extract_peak_Radon(GRad,Dbp,Tmax),
% figure
% imagesc(GRad_peaksc)
% pause
% 
% 
% figure
% imagesc(GRad_peaksr)
% pause

%%%%%%%%%%%%%%%%%
%l=L0;
for irho=1:s_rho
    %irho,
   for itheta=1:s_theta     
       %%%%%%%%%%%%%%%%%%%%%%%%%
       if(GRad_peaksc(irho,itheta) ~= 0)
           %ok=1,
           %pause
           %if(GRad(irho,itheta) ~= 0)
        startl=-nn/2;
        endl = -nn/2 + L0;
        while(startl <= nn/2)
            %while(startl <= nn/2 && endl <=nn/2)
            %ok=2,
            %pause
            for l=startl:endl
            %for l = 1:L0
            %ok=3,
            %pause
                [ixc] = find(xc == floor(rho(irho)*sintab(itheta)+l*costab(itheta))); %get indices of right xc values
          
                [iyr] = find(yr == floor(rho(irho)*costab(itheta)-l*sintab(itheta))); %get indices of right yr values
                %pause
                if(~isempty(ixc) && ~isempty(iyr))                    
                    Orig_ima(yr(iyr)+nn/2 + 1,xc(ixc)+nn/2 + 1) = 100;
                    %ok=4,
                    %pause
                end
            end
            startl = startl+L0;
            endl = endl + L0;
        end%while
        
       %else %make a straight segment 
           
       end%if 
       
       %%%%%%% smooth the current processed layer
%        if(GRad_peaksr(irho,itheta) ~= 0)
%            %ok=1,
%            %pause
%            %if(GRad(irho,itheta) ~= 0)
%         startl=-nn/2;
%         endl = -nn/2 + L0;
%         while(startl <= nn/2)
%             %while(startl <= nn/2 && endl <=nn/2)
%             %ok=2,
%             %pause
%             for l=startl:endl
%             %for l = 1:L0
%                 [ixc] = find(xc == floor(rho(irho)*sintab(itheta)+l*costab(itheta))); %get indices of right xc values
%           
%                 [iyr] = find(yr == floor(rho(irho)*costab(itheta)-l*sintab(itheta))); %get indices of right yr values
%                 %pause
%                 if(~isempty(ixc) && ~isempty(iyr))
%                     Orig_ima(yr(iyr)+nn/2+1,xc(ixc)+nn/2+1) = 100;
%                     %pause
%                 end
%             end
%             startl = startl+L0;
%             endl = endl + L0;
%         end%while
%         
%        %else %make a straight segment 
%            
%        end%if
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   end%itheta
end%irho
%%%%%%%%%%%%%%%%


% figure
% imagesc(Orig_ima)

end