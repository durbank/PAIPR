% Function to manually draw annual layers for PAIPR-processed echogram

function [man_layers] = draw_manual(radar)

% Preallocate cell array for position subscripts of manual layers
man_layers = cell(1,250);
i = 1;
draw = true;
while draw==true
    
    % Draw position of manual layers in radargram (layers should be continuous
    % across the entire radargram)
    f_draw = figure;
    imagesc(radar.data_smooth, [-3 3])
    hold on
%     for j = 1:length(guides)
%         plot(guides{j}(:,1), guides{j}(:,2), 'm')
%     end
    if i >= 2
        for k = 1:i-1
            plot(man_layers{k}(:,1), man_layers{k}(:,2), 'r')
        end
    end
    hi = drawpolyline();
    
    if isvalid(hi)
        
        % Find the range of the manually picked layer
%         col = (1:length(radar.Easting))';
        col = (max([1 round(min(hi.Position(:,1)))]):...
            min([round(max(hi.Position(:,1))) length(radar.Easting)]))';
        
        % Linearly interpolate layer row positions to the full range of the
        % manually picked layer
        row = interp1(hi.Position(:,1), hi.Position(:,2), col, ...
            'linear', 'extrap');
        
        % Export layer position subscripts to preallocated cell array
        man_layers{i} = [col row];
        i = i+1;
        
        close(f_draw)
        clear hi
    else
        draw = false;
    end
    
end

% Remove empty cells
man_layers = man_layers(~cellfun(@isempty,man_layers));

% Keep only traced layer data inside the boundaries of the echogram
keep_idx = cellfun(@(x) round(x(:,2))<=size(radar.data_smooth,1), ...
    man_layers, 'UniformOutput', false);
man_layers = cellfun(@(x,y) x(y,:), man_layers, keep_idx, ...
    'UniformOutput', false);


end