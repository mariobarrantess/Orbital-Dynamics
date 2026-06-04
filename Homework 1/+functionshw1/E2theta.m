function theta = E2theta (E, e)
% developed in computer class
    theta = acos((cos(E)-e)/(1-e*cos(E)));
end