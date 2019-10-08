classdef gen_node < handle
    
    
    properties  (SetAccess = private)%, Hidden = true)
        % cannot change:
        ind;
        r;
        phi;              
        m;
        strokeStrength;
        strokeDirection;  
        parent; 
        pAbs; % array 2x1 
        generatedForce; %array 2x1
        generatedTorque;
    end
    
%     properties         % CAN CHANGE
%     end
    
     methods
         function nd = gen_node(ind, node_mass, posXyz, radius, angle, parent)
             nd.ind = ind;
             nd.m = node_mass;
             nd.r = radius;
             nd.phi = angle;    
             nd.parent = parent;    
             nd.pAbs = posXyz;
             nd.strokeStrength = 0;
             nd.strokeDirection = 0; %0 is parallel to the body orientation
             nd.generatedForce = [0; 0; 0]; %vector 2x1
             nd.generatedTorque = [0; 0; 0];
         end 
         
         
         function setStrokeStrength(nd, str_stren)
             nd.strokeStrength = str_stren;
         end
         
         function setStrokeDirection(nd, str_dir)
             nd.strokeDirection = str_dir;
         end
         
         
         
         function p = get_abs_pos(nd)             
             p = [nd.parent.center(1) + nd.r * cosd(nd.phi + nd.parent.bodyAngle);...
                 nd.parent.center(2) + nd.r * sind(nd.phi +  nd.parent.bodyAngle); 0];
         end
         
         
         
         function  update(nd)
             % Recalculate node position and direciton of generated force.
             % Recalculate after each parent movement!
             
             % get abs position:
             nd.pAbs = nd.parent.facesAbs(nd.ind, :)';
           
             % STROKE:
             % recalculate direction of the stroke according to the new parent
             % position
             parentMainOrientaton       = nd.parent.bodyAngle;
             nodeRelativeOrientation    = nd.strokeDirection;
             nodeAbsoluteOrientation    = parentMainOrientaton + nodeRelativeOrientation;

             strokeStr                  = nd.strokeStrength;
             nd.generatedForce          = strokeStr*[ cosd(nodeAbsoluteOrientation); sind(nodeAbsoluteOrientation); 0];
             
             % TORQUE:
             force      = nd.generatedForce;
             radius     = nd.pAbs - nd.parent.center;
             T          = cross(radius, force);
             nd.generatedTorque = T;
         end
         
     end
end