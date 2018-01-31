function [out_layers] = smoothpicks_LTARE(ima_layers)
%% Smooth picked layers from LTARE automated pick
% this code will take the ima_layers output from the LTARE automated layer
% picker and will use a spline fit to smooth the layers and make them
% continuous over the entire echogram

%find the min and max layer numbers
layer_first = 1;
layer_last = max(max(ima_layers));
[rows,cols]=size(ima_layers);%v
out_layers = zeros(rows, cols);%v

% Fit the spline

for i=1:layer_last
    %i,
    
    [y_pick, x_pick]=find(ima_layers == i); %for each layer referenced by a number
    if length(x_pick) > size(ima_layers,2)/5 %want to eliminate any tiny apprearing line segments
              
       yy=interp1(x_pick,y_pick,1:cols,'nearest','extrap'); % interpolates to get smoothed y values
       
         for k=1:size(yy,2)
             if(fix(yy(k)) > 0 && fix(yy(k)) <= rows) %remember this an approximation, we may have negative values or values>rows
            out_layers((fix(yy(k))),k) =  i;             
             end
         end
    end
end


