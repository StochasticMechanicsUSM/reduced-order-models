function LocalAxes = BCLocalAxes(SpatialCoord,Nodes)
%% Beam/column local axes 
% Initial node
xi        = SpatialCoord(Nodes(1),:)';
% Final node
xj        = SpatialCoord(Nodes(2),:)'; 
% Local axes generation
% Axis 1
e1 = (xj-xi)/norm(xj-xi);
% Tolerance to define the axis 1 parallel to Z axis
Ztol = 10^-3;
if sqrt(e1(1)^2 + e1(2)^2) <= Ztol
    % Parallel to z
    % Axis 1
    e1 = [0;0;sign(e1(3))];
    % Axis 2
    e2 = [1;0;0];
    % Axis 3
    e3 = [0;sign(e1(3));0];
else
    % General case
    v3 = cross(e1,[0;0;1]);
    % Axis 3    
    e3 = v3/norm(v3);
    % Axis 2   
    e2 = cross(e3,e1);
end
% Local Axes (size -> 3x3)
LocalAxes = [e1';e2';e3'];
return
   