function runOnMovie(wrm, opts, movieStack)

textColor           = 'b';
bodyStampEachN      = 10;
init_pos            = wrm.center;
nFrames             = opts.simulation.frames;

xyzCollection       = wrm.center';
angleCollection     = wrm.bodyAngle;

timeFrameMat        = linspace(1, nFrames, size(movieStack,3));
img_W               = size(movieStack,2);

if opts.simulation.saveVideo    
    v_id = VideoWriter([opts.simulation.VideoName '.avi']); 
    % v_id.CompressionRatio = 3;
    open(v_id)
end



screenSize =  get(0,'screensize');
screenW = screenSize(3);
screenH = screenSize(4);

movie_frame1 = movieStack(:,:,1);
windowW = size(movie_frame1, 2);
windowH = size(movie_frame1, 1);

scaleWindow = min(screenH/windowH, screenW/windowW);
windowW = 0.8 * windowW * scaleWindow;
windowH = 0.8 * windowH * scaleWindow;

windowX = round((screenW - windowW)/2);
windowY = round((screenH - windowH)/2);

f1 = figure('position', [ windowX   windowY   windowW   windowH]);
colormap gray;
set(gca,'YDir','normal')
axis tight
hold on

  
%% deal with the body marks:
bodyMarkCounter     = 1;
bodyMarkCollection  = {};
bodyStampsColors    = jet(ceil(opts.simulation.frames / bodyStampEachN));

%% start
counter = 1;

while counter < nFrames
    
    counter = counter  + 1;
    if mod(counter, bodyStampEachN) == 0
        disp(counter)
    end
    
    wrm.move;    

    % collect data:
    angleCollection     = [angleCollection; wrm.bodyAngle];
    xyzCollection       = [xyzCollection; wrm.center'];
    
    [~,findMovieFrame] = min(abs(timeFrameMat - counter));
    if findMovieFrame>size(movieStack,3)
        findMovieFrame=size(movieStack,3);
    end
        
    figure(f1);
    cla;
    xlim([1, img_W])
    imagesc(movieStack(:, :, findMovieFrame))
    plot(init_pos(1), init_pos(2), 'rx', 'LineWidth', 2)
    plot(xyzCollection(:, 1), xyzCollection(:, 2), 'm-') 

    wrm.addToFig(f1);
    
    centrCoord = [wrm.center(1), wrm.center(2)];
    %text(10,10, ['frame  ' num2str(counter)], 'color', textColor)
    text(10,10, ['iteration  ' num2str(counter)], 'color', textColor)
    text(10,20, ['stroke  ' num2str(wrm.strokeForce)], 'color', textColor)
    text(10,30, ['angle  ' num2str(wrm.bodyAngle)], 'color', textColor)
    text(10,40, ['center  ' num2str(centrCoord)], 'color', textColor)
    
    drawnow

    
    
        if opts.simulation.saveVideo
            frame_i = getframe(f1);
            writeVideo(v_id, frame_i)   
        end
 
       
    
end
close(v_id)
