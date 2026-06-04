% function theta = E2theta (E, e)
%     theta = acos((cos(E)-e)/(1-e*cos(E)));
% end


function theta = E2theta(E, e)
theta = atan2( sqrt(1-e^2)*sin(E), cos(E)-e );
if theta < 0
    theta = theta + 2*pi;
end
end