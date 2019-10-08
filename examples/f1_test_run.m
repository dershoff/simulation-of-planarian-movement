function  f1_test_run

% INITIALIZE A WORM:
opts.worm.mutant = 'vfl3'; 
% 'vfl3' - left
% 'odf2' - right
% 'wt'
opts.worm.length = 8; %mm
opts.worm.bodyRatio = 2.8; %length to width, normal value = 2.8

opts.worm.initialPosition = [0, 0, 0];
opts.worm.orientation = 90; %90 = upwards
opts.worm.strokeForce = 0.2; %arbitrary unit
opts.worm.mass = 10; %arbitrary
opts.worm.dt = 1; %time step for integration
w = gen_worm(opts);

opts.simulation.frames = 100;
opts.simulation.showFrames = 1;
opts.simulation.showStamps = 0;
opts.simulation.saveVideo = 0;
opts.simulation.VideoName = 'temp_test';

plotOn =  1;
runStuff(w, opts, plotOn);
