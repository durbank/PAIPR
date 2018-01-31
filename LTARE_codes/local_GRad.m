function [GRad]=local_GRad(ima_in,L0,theta)
%calling: [GRad]=local_GRad(ima_in,L0)
%input parameters: ima_in: input imaged echogram
%                  L0    : intergration length, typically in [1, 5] pixels
%                  theta : angular direction range typicall theta in [-pi/5,
%                  pi/5] for ku-band radar and in [-pi/64, pi/64] for radar
%                  depth sounder
%
%output parameters: Grad: result radon transform domain
%                   
%
%Description: local_GRad computes the forward radon transform of ima_in 
%using L0 integration length, and the angular direction range theta.


%%%%%%%%%% Initialization %%%%%%%%%%%%%
[nr,nc] = size(ima_in);
nn=min(nr,nc);
%pause
ynr = floor((nn+1)/2);
xnc = floor((nn+1)/2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = -nn/2:nn/2-1;
s_rho = size(rho,2); % nn points

%theta = -pi/6:pi/(3*nn):pi/6-pi/(3*nn); %ok 
%theta = -pi/5:pi/(2.5*nn):pi/5-pi/(2.5*nn); %ok+ 
%theta = -pi/4:pi/(2*nn):pi/4-pi/(2*nn); 
%theta = -pi/3:pi/(1.5*nn):pi/3-pi/(1.5*nn);  
%theta = -pi/2:pi/nn:pi/2-pi/nn; 

costab = cos(theta);
sintab = sin(theta);
s_theta = size(theta,2); %nn points

GRad = zeros(s_rho,s_theta);

%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%

for irho=1:s_rho
    %irho,
   for itheta=1:s_theta
              
       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%        mmin = min(rho(irho)*costab(itheta) - sig*sintab(itheta), rho(irho)*costab(itheta)-(sig+l)*sintab(itheta)); 
%        mmin = ceil(mmin)+nn/2;
%        mmax = max(rho(irho)*costab(itheta) - sig*sintab(itheta), rho(irho)*costab(itheta)-(sig+l)*sintab(itheta));
%        mmax = ceil(mmax)+nn/2;
%        
%        nmin = rho(irho)*sintab(itheta) + sig*costab(itheta); %==y
%        nmin = ceil(nmin)+nn/2;
%        nmax = rho(irho)*sintab(itheta) + (sig+l)*costab(itheta);
%        nmax = ceil(nmax)+nn/2;
%              
%        
%        for m=mmin:mmax
%            for n=nmin:nmax
%                if(rho(irho) == floor(m*costab(itheta) + n*sintab(itheta)))
%                 if(n>0 && n<=nn && m>0 && m<=nn)                                
%                     GRad(irho,itheta) = GRad(irho,itheta) + ima_in(m,n);                    
%                 end
%                end
%            end
%        end
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      S_L0 = 0; %first sum 
      startl=-nn/2;
      endl = -nn/2 + L0;
      for l=startl:endl
         xc = rho(irho)*sin(theta(itheta)) + l*cos(theta(itheta)) + nn/2; xc = floor(xc);
         yr = rho(irho)*cos(theta(itheta)) - l*sin(theta(itheta)) + nn/2; yr = floor(yr);
         if(xc>0 && xc<=nn && yr>0 && yr<=nn)                                
            %GRad(irho,itheta) = GRad(irho,itheta) + ima_in(yr,xc); 
            S_L0 = S_L0 + ima_in(yr,xc); 
         end
      end%l
      
      GRad(irho,itheta) = GRad(irho,itheta) + S_L0; 
      while(endl <= nn/2)
         % while(startl < nn/2 && endl < nn/2)
        l=startl;    
        xc = rho(irho)*sin(theta(itheta)) + l*cos(theta(itheta)) + nn/2; xc = floor(xc);
        yr = rho(irho)*cos(theta(itheta)) - l*sin(theta(itheta)) + nn/2; yr = floor(yr);
        if(xc>0 && xc<=nn && yr>0 && yr<=nn)
            S_L0 = S_L0 - ima_in(yr,xc); 
        end
        l = endl + 1;
      
        xc = rho(irho)*sin(theta(itheta)) + l*cos(theta(itheta)) + nn/2; xc = floor(xc);
        yr = rho(irho)*cos(theta(itheta)) - l*sin(theta(itheta)) + nn/2; yr = floor(yr);
        if(xc>0 && xc<=nn && yr>0 && yr<=nn)
            S_L0 = S_L0 + ima_in(yr,xc); 
        end
      
        GRad(irho,itheta) = GRad(irho,itheta) + S_L0; %update summation along the line
        startl = startl + 1;
        endl = endl + 1;
      end%while 
      GRad(irho,itheta) = GRad(irho,itheta);%/nn;
      
   end%itheta
end%irho
% figure
% imagesc(ima_in(1:nn,1:nn))
% figure
% imagesc(GRad)

end