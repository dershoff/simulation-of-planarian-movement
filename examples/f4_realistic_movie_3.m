function  f4_realistic_movie_3

movie_path      = 'D:\worms\vfl3\';
movieName       = 'Vfl3';
startFrame      = 135;
endFrame        = 233;
framesPerSecond = 351/53;               %total frames/total sec
timeStep        = 1/framesPerSecond;    % sec
window_size_pix = [644, 358];           % as in imageJ
window_size_um  = [21207.23, 11789.11];   % as in imageJ
simTimeStep     = 0.01;                 % 10 ms

% measured parameters:
x = 5729.88;        %um
y = 4149.35;        %um
L = 10000;          %um
w = 1200.95;         %um
a = 32;        %deg

pix_to_um = mean(window_size_um./window_size_pix);
um_to_pix = 1/pix_to_um;

% INITIALIZE A WORM:
% 'vfl3' - left
% 'odf2' - right
% 'wt'
opts.worm.mutant            = 'vfl3'; 


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
opts.simulation.VideoName   = 'WT_test_2';

wrm = gen_worm(opts);
%wrm.singleShow;




movieStack = [];
for i = startFrame : endFrame
    movieStack = cat(3, movieStack, imread([movie_path, movieName, sprintf('%0.4d', i) ,'.tif']));
end

close all;
runOnMovie(wrm, opts, movieStack);



   

