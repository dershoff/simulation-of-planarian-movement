function [xyzCollection, angleCollection] = runStuff(wrm, opts, plotResultON)

% input:
% w - is instance of the worm class;

bodyStampEachN      = 50;
init_pos            = wrm.center;
nFrames             = opts.simulation.frames;

xyzCollection       = wrm.center';
angleCollection     = wrm.bodyAngle;


if opts.simulation.showFrames
    
    close all    
    f1 = figure('position', [      142         299        1446         608]);
    subplot(1,2,1)
    xlim([-200 200])
    ylim([-200 200])
    % axis equal
    hold on;
    
    subplot(1,2,2)
    % view(30,25)
    view(2)
    axis equal
    hold on;
end

if opts.simulation.saveVideo
    
    v_id = VideoWriter([opts.simulation.VideoName '.avi']); 
    % v_id.CompressionRatio = 3;
    open(v_id)
end



%% deal with the body marks:
bodyMarkCounter     = 1;
bodyMarkCollection  = {};
bodyStampsColors    = jet(ceil(opts.simulation.frames / bodyStampEachN));

%% start
done = 0;
counter = 0;

while ~done
    
    wrm.move;    
    counter = counter  +1;

    % collect data:
    angleCollection     = [angleCollection, wrm.bodyAngle];
    xyzCollection       = [xyzCollection; wrm.center'];
    
    if mod(counter, bodyStampEachN) == 0
        disp(counter)
    end
     
     
    if opts.simulation.showFrames

        subplot(1,2,1)
        cla;
        plot(init_pos(1), init_pos(2),'rx', 'LineWidth', 2) 
        plot(xyzCollection(:, 1), xyzCollection(:, 2), 'k.') 
        plot(wrm.center(1), wrm.center(2), 'r.')

        if opts.simulation.showStamps
            if mod(counter, bodyStampEachN) == 0
                col                 = bodyStampsColors(bodyMarkCounter, :);
                bodyMark            = wrm.stamp(col);
                bodyMarkCollection  = [bodyMarkCollection, {bodyMark}];
                bodyMarkCounter     = bodyMarkCounter + 1;
            end
            for mmm = 1 : length(bodyMarkCollection)
                bm = bodyMarkCollection{mmm};
                plot(bm(:,4), bm(:,5), 'r.')
                plot(bm(:,1), bm(:,2), '-', 'color', bm(1, 7:9), 'LineWidth', 2)
            end
        end
              
        subplot(1,2,2)
        cla;
        wrm.show;    
        title(['Frame : ' num2str(counter) '; Body angle: ' num2str(wrm.bodyAngle) '; Torque: ' num2str(norm(wrm.totalTorque)) ])
        % title(['Body angle: ' num2str(w.oAngle) '; Stroke Dir: ' num2str(StrokeAngleDefault)])
            
        drawnow;
    end     

    
    
    if counter <= nFrames
        
        if opts.simulation.saveVideo
            frame_i = getframe(f1);
            writeVideo(v_id, frame_i)           
        end
    else
         done = 1;
         if opts.simulation.saveVideo
             close(v_id)
         end        
    end
    
end



if plotResultON
%close all

figure('position', [      142         299        1446         608]);
subplot(1,2,1)
title(wrm.mutant)

axis equal
hold on;
plot(init_pos(1),init_pos(2), 'rx', 'LineWidth', 2)
plot(xyzCollection(:,1), xyzCollection(:,2),'k-') 

xl = xlim;
yl = ylim;
xlim(1.2*xl);
ylim(1.2*yl);
    
subplot(1,2,2)
title('Angle evolution. 90 = upwards.')
xlabel('time (iterations)')
ylabel('angle (degree)')
hold on
plot(angleCollection)
ylim([45 135])
end