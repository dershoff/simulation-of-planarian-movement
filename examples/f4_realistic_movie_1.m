function  f4_realistic_movie_1

movie_path      = 'D:\worms\wt\1\';
movieName       = 'WT-10fps';
startFrame      = 88;
endFrame        = 146;
framesPerSecond = 250/21;               %total frames/total sec
timeStep        = 1/framesPerSecond;    % sec
window_size_pix = [518, 258];           % as in imageJ
window_size_um  = [11371.87, 5663.98];  % as in imageJ
simTimeStep     = 0.01;                 % 10 ms

% measured parameters:
x = 2480.74;        %um
y = 2831.9;         %um
L = 4282.210;       %um
w = 740;            %um
a = -14;          %deg

pix_to_um = mean(window_size_um./window_size_pix);
um_to_pix = 1/pix_to_um;

% INITIALIZE A WORM:
% 'vfl3' - left
% 'odf2' - right
% 'wt'
opts.worm.mutant            = 'wt'; 


opts.worm.length            = L * um_to_pix;    % pix
opts.worm.bodyRatio         = L/w;              % length to width, normal value = 2.8
opts.worm.initialPosition   = [x, y, 0] * um_to_pix;
opts.worm.orientation       = a;                % from ImageJ
opts.worm.strokeForce       = 0.032;              % arbitrary unit
opts.worm.mass              = 1;               % arbitrary
opts.worm.dt                = 1;                % time step for integration

opts.simulation.frames      = ceil((endFrame - startFrame)*timeStep/simTimeStep);       % number of frames
opts.simulation.showFrames  = 1;
opts.simulation.showStamps  = 1;
opts.simulation.saveVideo   = 1;
opts.simulation.VideoName   = 'WT_test_1';

wrm = gen_worm(opts);
%wrm.singleShow;




movieStack = [];
for i = startFrame : endFrame
    movieStack = cat(3, movieStack, imread([movie_path, movieName, sprintf('%0.4d', i) ,'.tif']));
end

% movie_frame1 = movieStack(:,:,1);
% windowW = 2*size(movie_frame1, 2);
% windowH = 2*size(movie_frame1, 1);
% 
% f1 = figure('position', [ 333   326   windowW   windowH]);
% colormap gray;
% set(gca,'YDir','reverse')
% axis tight
% hold on
% 
% imagesc(movie_frame1)
% 
% wrm.addToFig(f1);



close all;

runOnMovie(wrm, opts, movieStack);



   

