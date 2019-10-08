function setNewAnglesForces(w, forceVal, strokeRelAngle)



w.singleShow;


[xLims, yLims] = ginput(2);
indList = [];
for n = w.nodeList
    x = n.pAbs(1); y = n.pAbs(2); z = n.pAbs(3);
    if (x>=min(xLims) && x<=max(xLims)) && (y>=min(yLims) && y<=max(yLims))
        indList = [indList, n.ind];
        plot3(x, y, z, 'co')
    end
end
% StrokeAngleToSet = input('enter the angle for the selected nodes: ');
% StrokeForceToSet = input('enter the force for the selected nodes: ');

StrokeAngleToSet = strokeRelAngle;
StrokeForceToSet = forceVal;

StrokeForceList = StrokeForceToSet*ones(1, length(indList));
StrokeAngleList = StrokeAngleToSet * ones(1, length(indList));
w.setNodeStrokes(indList, StrokeForceList, StrokeAngleList)

% show initial state:
w.singleShow;
