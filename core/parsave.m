% A function to save the output of radar processing within a parellel for
% loop (so as to preserve variable transparency)

function [success] = parsave(mdata, output)

save(output, '-struct', 'mdata', '-v7.3')

success = true;

end