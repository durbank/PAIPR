function [nppy]=find_last_y_nz(ppy)
%the function finds the last y values within breaks 
% and return the corresponing new ppy called nppy where breaks are filled
% with the mean of previous y values


[iy]=find(ppy == 0);

nppy = ppy;
if(~isempty(iy))
    for j=1:length(iy)
     i=iy(j);
     while(i <= length(ppy) && nppy(i) == 0)
       nppy(i) = mean(ppy(1:iy(j)-1));
       i=i+1;
     end
    end
end
end