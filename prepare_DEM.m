function prepare_DEM(DEM_fname, SHAPE_OUT_FNAME, TOPOTOOLBOX_PATH, Area_threshold)
%generates flow

% Example parameters for test run:
%TOPOTOOLBOX_PATH = '/home/bodo/Dropbox/Matlab-work/topotoolbox';
%DEM_fname='/home/bodo/Dropbox/Himalaya/ChanGeom_Himalaya/WHimalaya/WHimalaya_srtm1_UTM44N_WGS84.tif';
%shape_fname='/home/bodo/Dropbox/Himalaya/ChanGeom_Himalaya/WHimalaya/WHimalaya_srtm1_UTM44N_WGS84_STR_1e4m.shp';
%Area_threshold = 1e4;
[pathstr,fname,ext] = fileparts(DEM_fname);
%
%make sure topotoolbox is in the PATH, e.g.:
addpath(genpath(TOPOTOOLBOX_PATH));
%

DEM = GRIDobj(DEM_fname);
FD  = FLOWobj(DEM,'preprocess','carve','mex',true);
DEM = imposemin(FD,DEM,0.0001);
A   = flowacc(FD);
area_1km_in_pixels = round(Area_threshold/(A.cellsize^2));
w = A > area_1km_in_pixels;
S = STREAMobj(FD,w);
D_maxfromchannel = distance(S, 'max_from_ch');
g = gradient(S,DEM);
AOI_DEM_gradient8 = gradient8(DEM);
% calculate ksn with theta = 0.45
theta = -0.45;
AOI_ksn045 = A;
AOI_ksn045 = AOI_DEM_gradient8 ./ (A.*...
    (A.cellsize^2)).^theta;
AOI_ksn045.name = 'ksn045';
AOI_ksn045.zunit = 'm^0.9';

fig_fname = sprintf('%s%s_DEM.jpg', pathstr, fname);
if exist(fig_fname,'file') ~= 2
    figure('visible','off'), clf
    fprintf(1,'generating overview figure\n');
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', 'A4');
    imageschs(DEM, DEM, 'ticklabels', 'nice')
    grid on
    hold;
    plot(S, 'color', [0.5 0.5 0.5], 'Linewidth', 0.5)
    s_string = sprintf('DEM: %s', fname);
    title(s_string, 'Fontsize', 16)
    if exist('export_fig') == 2
        export_fig(fig_fname)
    else
        saveas(gcf,fig_fname, 'jpg');
    end
end

fig_fname = sprintf('%s%s_ksn.jpg', pathstr, fname);
if exist(fig_fname,'file') ~= 2
    figure('visible','off'), clf
    fprintf(1,'generating ksn figure\n');
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', 'A4');
    imageschs(DEM, AOI_ksn045, 'ticklabels', 'nice', 'caxis', [0 200])
    hold;
    s_string = sprintf('K_{sn}: %s', fname);
    title(s_string, 'Fontsize', 16)
    if exist('export_fig') == 2
        export_fig(fig_fname)
    else
        saveas(gcf,fig_fname, 'jpg');
    end
end

fprintf(1,'Writing to shapefile: %s\n', SHAPE_OUT_FNAME);
shapewrite(STREAMobj2mapstruct(S, 'attributes',{'elev' DEM @mean 'A_km2' A*A.cellsize^2/1e6 @mean 'flowacc' A @mean 'gradient' g @mean 'd_maxdd' D_maxfromchannel @mean 'd_out' S.distance @mean 'ksn' AOI_ksn045 @mean}), SHAPE_OUT_FNAME);
exit
