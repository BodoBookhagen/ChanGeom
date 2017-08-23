function chanextract(inputtif, exporttif, start_pt) 
% function chanextract(inputtif, exporttif, start_pt) 
%
% Input parameters
% <inputtif>:   logical or binary TIF file. 
% <exporttif>:  name of the exported channel width data geotiff. This file will have width values 
% <cellsize>:   taken from <inputtif> TIF file
% <start_pt>:   0 for processing from north to south (i.e., start point is to the south of the end point)
%               1 for processing from south to north (i.e., start point is to the north of the end point)
%
% Written by Burch Fisher and Bodo Bookhagen, May 2014
%
% When using this code, cite: 
% Fisher, G.B., Bookhagen, B., Amos, C.B. (2013): "Channel planform geometry
% and slopes from freely available high-spatial resolution imagery and DEM
% fusion: Implications for channel width scalings, erosion proxies, and
% fluvial signatures in tectonically active landscapes", Geomorphology,
% doi:10.1016/j.geomorph.2013.04.011    
%
%

%by default, Matlab save centerlines as mat files in the following
%directory:
output_mat_dir = 'MAT_centerline';
%find directory names
[foo1, foo2] = strtok(exporttif,'.');
if isunix() == 1 || ismac() == 1
    slashc = '/';
    slashregexp = '/*.tif';
elseif ispc() == 1
    slashc = '\';
    slashregexp = '/*.tif';
end
foo1_id = regexp(foo1,slashc);
if length(foo1_id) > 0
    foo_dir = exporttif(1:foo1_id(end-1));
    slash_id = regexp(exporttif, slashc);
    tif_id = regexp(exporttif, slashregexp);
    matfile2save = exporttif(slash_id(end)+1:tif_id-1);
else
    %no / in dir name
    foo_dir='';
    matfile2save = strtok(exporttif,'.');
end
MATname = strcat(foo_dir, output_mat_dir,slashc, matfile2save);
full_output_mat_dir=strcat(foo_dir, output_mat_dir);
clear slash_id tif_id foo1 foo2 foo_dir foo1_id
if exist(full_output_mat_dir) == 0
    mkdir(full_output_mat_dir)
end
if exist(inputtif) == 0
    fprintf('%s: file not found\n', inputtif);
    exit;
end

[riv_msk_fil,header] = geotiffread(inputtif);
cellsize = header.DeltaX;
info = geotiffinfo(inputtif);
if info.BitDepth ~= 1
    fprintf('%s: not in logical (or 1-bit) format\n', inputtif);
    exit;
end
%Extract NoData value, usually the mode of the data (because most of the
%data are background data
%Alternatively, could use gdalinfo and look for NoDATA entry.
NoData = mode(single(riv_msk_fil(:)));

%Need to make outline polygon from tiff file
riv_msk = bwmorph(riv_msk_fil, 'remove');
[idx1, idx0, nr1, riv_msk_fil] = riv_mask_fil(riv_msk);

%Create centerline
fprintf(1,'%s: filtering binary grid...\n', inputtif);
centerline = bwmorph(riv_msk_fil, 'thin', 500);
centerline = bwmorph(centerline, 'clean', 10);

%Prune with endpoints until only two endpoints exist (start & beginning)
no_end_points = centerline;nr_of_endpoints = 3; nr_of_iterations = 0;
while(nr_of_endpoints > 2)
    nr_of_iterations = nr_of_iterations + 1; end_points = endpoints(no_end_points); no_end_points = no_end_points & ~end_points; nr_of_endpoints = length(find(end_points > 0));
    if nr_of_iterations > 10000, fprintf(1,'%s: more than 10,000 iterations in finding endpoints - returning\n', inputtif); return; end
end
centerline = no_end_points; clear no_end_points clear nr_of_endpoints nr_of_iterations

%Alternative approach to prune spurs, but slow:
%centerline = bwmorph(centerline, 'spur', 100); %this takes a long time and may have to be altered - it removes start and endpoints, too

centerline = bwmorph(centerline, 'thin');
idx_centerline = find(centerline == 1); idx_no_centerline = find(centerline == 0);

%find southernmost point
centerline_xy = size(centerline); 
end_points = endpoints(centerline); [endpoint_idxx, endpoint_idxy] = find(end_points > 0);

%first coordinate in centerline_xy is Y extend - find coordinate in
%endpoint_idxx that is closest to centerline_xy(1)
idx_south = find(min(centerline_xy(1) - endpoint_idxx) == (centerline_xy(1) - endpoint_idxx));
idx_north = find(max(centerline_xy(1) - endpoint_idxx) == (centerline_xy(1) - endpoint_idxx));
south_point = [endpoint_idxx(idx_south) endpoint_idxy(idx_south)];
north_point = [endpoint_idxx(idx_north) endpoint_idxy(idx_north)];
centerline_idx = bwtraceboundary(centerline, south_point, 'N', 8);
if isempty(centerline_idx), fprintf(1,'%s: Can not find bwtraceboundary of centerline\n', inputtif); return; end

%centerline_idx no contains the pixel coordinates of the filtered centerline
clear idx_south idx_north end_points endpoint_idxx endpoint_idxy

%Using 'Quasi-Euclidean Distance' to calculate distance from edge to centerline and multiply by 2    
riv_msk_distance_qe = single(bwdist(riv_msk, 'quasi-euclidean')); riv_msk_distance_qe(idx0) = 0; riv_msk_distance_qe_idx = find(riv_msk_distance_qe == 1);
rwidth1 = riv_msk_distance_qe; rwidth1(idx_no_centerline) = 0; rwidth1 = single(rwidth1 .* (2 * cellsize)+cellsize); rwidth1(idx_no_centerline) = 0;%clear riv_msk_distance
%This adds one pixel to the width and then resets the background values to 0
%save riverwidth values
rwidth1_col = single(NaN(size(centerline_idx,1),1));
for i = 1:size(centerline_idx,1), rwidth1_col(i) = rwidth1(centerline_idx(i,1), centerline_idx(i,2)); end
%rwdith1 contains riverwidth in meters along the centerline from south to north

index0 = find(rwidth1==0);
rwidth1(index0)= NoData;

%Save Geotiff of Rwidth Result
key=info.GeoTIFFTags.GeoKeyDirectoryTag;
geotiffwrite(exporttif,rwidth1, header, 'GeoKeyDirectoryTag', key);

%% In this section, you will calculate the cumulative distance along the channel
% Find the centerline distances
if start_pt == 0; start = south_point;
else start = north_point;
end    

% Define first_pixelx and first_pixely, the algorithm goes in any direction
% but you must have the coordinates of the last or first point along the
% centerline
first_pixely= start(1,2); %actually x value
first_pixelx = start(1,1); %actually y value 

centerline_numbers = single(NaN(size(centerline)));
centerline_distance = single(NaN(size(centerline)));
centerline_cumdist = single(NaN(size(centerline)));
centerline_processed = zeros(size(centerline)); centerline_processed = logical(centerline_processed);

runningdistance = 0;
counter = 1;
whilecounter = 1;
stop = 0;

while (stop == 0)
    if (mod(whilecounter, 1000) == 0), fprintf('%d, ', whilecounter); end
    if (mod(whilecounter, 20000) == 0), fprintf('\n'); end
    % start at first_pixelx, first_pixely
    if whilecounter == 1
        centerline_numbers(first_pixelx, first_pixely) = counter; centerline_distance(first_pixelx, first_pixely) = 0;
        centerline_cumdist(first_pixelx, first_pixely) = 0; centerline_processed(first_pixelx, first_pixely) = 1;
        if ((centerline(first_pixelx,first_pixely-1) == 1) & (centerline_processed(first_pixelx,first_pixely-1) == 0)), next_x = first_pixelx; next_y = first_pixely-1; step_y = 1; step_x = 0;
        elseif ((centerline(first_pixelx,first_pixely+1) == 1) & (centerline_processed(first_pixelx,first_pixely+1) == 0)), next_x = first_pixelx; next_y = first_pixely+1; step_y = 1; step_x = 0;
        elseif ((centerline(first_pixelx-1,first_pixely-1) == 1) & (centerline_processed(first_pixelx-1,first_pixely-1) == 0)), next_x = first_pixelx-1; next_y = first_pixely-1; step_y = 1; step_x = 1;
        elseif ((centerline(first_pixelx-1,first_pixely) == 1) & (centerline_processed(first_pixelx-1,first_pixely) == 0)), next_x = first_pixelx-1; next_y = first_pixely; step_y = 0; step_x = 1;
        elseif ((centerline(first_pixelx-1,first_pixely+1) == 1) & (centerline_processed(first_pixelx-1,first_pixely+1) == 0)), next_x = first_pixelx-1; next_y = first_pixely+1; step_y = 1; step_x = 1;
        elseif ((centerline(first_pixelx+1,first_pixely-1) == 1) & (centerline_processed(first_pixelx+1,first_pixely-1) == 0)), next_x = first_pixelx+1; next_y = first_pixely-1; step_y = 1; step_x = 1;
        elseif ((centerline(first_pixelx+1,first_pixely) == 1) & (centerline_processed(first_pixelx+1,first_pixely) == 0)), next_x = first_pixelx+1; next_y = first_pixely; step_y = 0; step_x = 1;
        elseif ((centerline(first_pixelx+1,first_pixely+1) == 1) & (centerline_processed(first_pixelx+1,first_pixely+1) == 0)), next_x = first_pixelx+1; next_y = first_pixely+1; step_y = 1; step_x = 1;
        else stop = 1; %no adjacent pixel found
        end
        distance = sqrt(step_x^2 + step_y^2);
        runningdistance = runningdistance + distance;
        centerline_distance(next_x,next_y) = distance;
        centerline_cumdist(next_x,next_y) = runningdistance;
        counter = counter + 1;
        centerline_processed(next_x, next_y) = 1;
    else
        %now find adjacent pixel
        if ((centerline(next_x,next_y-1) == 1) & (centerline_processed(next_x,next_y-1) == 0)), next_x = next_x; next_y = next_y-1; step_y = 1; step_x = 0;
        elseif ((centerline(next_x,next_y+1) == 1) & (centerline_processed(next_x,next_y+1) == 0)), next_x = next_x; next_y = next_y+1; step_y = 1; step_x = 0;
        elseif ((centerline(next_x-1,next_y-1) == 1) & (centerline_processed(next_x-1,next_y-1) == 0)), next_x = next_x-1; next_y = next_y-1; step_y = 1; step_x = 1;
        elseif ((centerline(next_x-1,next_y) == 1) & (centerline_processed(next_x-1,next_y) == 0)), next_x = next_x-1; next_y = next_y; step_y = 0; step_x = 1;
        elseif ((centerline(next_x-1,next_y+1) == 1) & (centerline_processed(next_x-1,next_y+1) == 0)), next_x = next_x-1; next_y = next_y+1; step_y = 1; step_x = 1;
        elseif ((centerline(next_x+1,next_y-1) == 1) & (centerline_processed(next_x+1,next_y-1) == 0)), next_x = next_x+1; next_y = next_y-1; step_y = 1; step_x = 1;
        elseif ((centerline(next_x+1,next_y) == 1) & (centerline_processed(next_x+1,next_y) == 0)), next_x = next_x+1; next_y = next_y; step_y = 0; step_x = 1;
        elseif ((centerline(next_x+1,next_y+1) == 1) & (centerline_processed(next_x+1,next_y+1) == 0)), next_x = next_x+1; next_y = next_y+1; step_y = 1; step_x = 1;
        else stop = 1; %no adjacent pixel found
        end
        
        %calculate distance
        distance = sqrt(step_x^2 + step_y^2);
        runningdistance = runningdistance + distance;
        centerline_distance(next_x,next_y) = distance;
        centerline_cumdist(next_x,next_y) = runningdistance;
        centerline_processed(next_x, next_y) = 1;
        centerline_numbers(next_x, next_y) = counter;
        counter = counter + 1;
    end
    whilecounter = whilecounter + 1;
end % while
fprintf('\n')

% Multiplies results times cellsize
centerline_distance = centerline_distance.*cellsize;
centerline_cumdist = centerline_cumdist.*cellsize;

% Converts NaN values in background to whatever the NoData value is
% specified as in the function
poop = isnan(centerline_cumdist);
centerline_cumdist(poop)= NoData;
centerline_numbers(poop)= NoData;

% Images the centerline_cumdist to make sure it didn't get truncated
% imagesc(centerline_cumdist); colorbar

% Uses the name of the .mat file and concatenates with the numbers and
% cumdist files for export
name = strtok(exporttif,'.');
cent_cumdist_name = strcat(name,'_cumdist.tif');

% EXPORT Geotiffs of Results
geotiffwrite(cent_cumdist_name,centerline_cumdist, header, 'GeoKeyDirectoryTag', key);
clear centerline_idx centerline_xy i idx0 idx1 idx_centerline idx_no_centerline index0 nr1...
    riv_msk_distance_qe_idx rwidth1_col cellsize riv_msk_cellsize full_output_mat_dir


clear counter distance next_x next_y poop runningdistance step_x step_y stop...
    whilecounter x y centerline_distance first_pixelx first_pixely centerline_processed...
    chanextractMATfile exporttif name cent_cumdist_name start key

save(MATname, '-v7.3'); 
exit;
