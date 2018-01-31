function [fppy,fppx] = clean_zeros_py(py,px)
%calling: [fppy,fppx] = clean_zeros_py(py,px)
%input parameters:  py: y values from a layer
%                   px: x values from a layer
%
%output parameters: fppy: py smoothed
%                   fppx: px smoothed 
%Description: This function clean zeros on layer coordinates (py,px)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ry = find(py); %get non zeros row indexes
rx = find(px); %get non zeros row indexes

fppy=[]; fppx=[];
ppy=[]; ppx=[];
%sry=size(ry),
if(~isempty(ry) && length(ry) > 1)      
    ppx = px;%(1,ry); %ok get valid px values
    [nppy]=find_last_y_nz(py);
    %ppy = py(1,ry); %ok get valid py values              
    ppy = nppy;  

    
    %ima_layers = zeros(400,1050);
    %ppx,
    %ppy,
 %==================================================================
%         y_pick_spline=spline(ppx,ppy); % the return value is a structure not a scalar nor a vector
%         yy=ppval(y_pick_spline,1:size(ppx,2)); %get the smoothed y values
%         fppy = fix(yy);
%         fppx= fix(ppx);
 %=====================================================================
  
%          for k=1:size(yy,2)
%              if(fix(yy(k)) > 0 && fix(yy(k)) <= size(ima_layers,1)) %remember this an approximation, we may have negative values or values>rows
%                 ima_layers((fix(yy(k))),k) =  1;             
%              end
%          end
      
%          figure
%         imagesc(ima_layers);
%          pause
       
    
 %+++++++++++++++++++++++++++better +++++++++++++++++++++++++++++++++++++++   
% ppx = smooth(ppx,'rlowess'),
% ppy = smooth(ppy,'rlowess'),

%smooth1
sx=1:1:length(ppx);
sy=ppy; sy=round(sy);
%zsy = ku_smoothn(abs(ppy),abs(ppx),'robust');
zsy = ku_smoothn(abs(sy),abs(sx),'robust');
%zsy = fix(zsy);
%sx = abs(ppx);
      zsy = ceil(zsy);
      %sx = ceil(sx);      
            
      fppx = sx;
      fppy = zsy;


  end
end