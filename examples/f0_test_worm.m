%% INITIALIZE A WORM:
opts.simulation.frames = 1000;
opts.simulation.showFrames = 0;
opts.simulation.saveVideo = 0;

opts.worm.length = 8; %mm
opts.worm.bodyRatio = 2.8; %length to width, normal value = 2.8
opts.worm.initialPosition = [5, 3, 0];
opts.worm.orientation = 90; %90 = upwards

opts.worm.strokeForce = 0.2; %arbitrary unit
opts.worm.mass = 10; %arbitrary
opts.worm.dt = 1; %time step for integration

opts.worm.mutant = 'vfl3'; 
% 'vfl3' - left
% 'odf2' - right
% 'wt'

% w = genWorm2(opts);
w = gen_worm(opts);

w.singleShow;