function theta = E2theta (E, e)
    theta = acos((cos(E)-e)/(1-e*cos(E)));
end