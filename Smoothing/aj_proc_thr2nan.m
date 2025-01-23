function fn_out = aj_proc_thr2nan(fn_in, threshold)
% Toy function to replace all values below a threshold contained in fn_in
% to NaN.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Prepare input/output
% Images are overwritten so same output file names
fn_out = fn_in;

% Turn a cell array into char array
if iscell(fn_in)
    fn_in = char(fn_in);
end

% Number of images to deal with
Nimg = size(fn_in,1);

%% Deal with all the images.
% Since current computer have large memories, we can proceed directly by
% using the file_array object to replace threshold>'s with NaN's, when needed.
for ii=1:Nimg
    V_ii = spm_vol(deblank(fn_in(ii,:))); % memory map each image
    if spm_type(V_ii.dt(1),'nanrep')
        % Do the conversion if it has NaN representation
        dd = V_ii.private.dat(:,:,:); % read in the whole 3D image
        dd(dd(:)<=threshold) = NaN; % Turn threshold>'s into NaN in data array
        V_ii.private.dat(:,:,:) = dd; % write back data into file
    end
end

end