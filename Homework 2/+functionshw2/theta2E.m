function E = theta2E (theta, e)
    E = acos((cos(theta)+e)/(1+e*cos(theta)));

    if mod(theta, 2*pi)>pi
        E = 2*pi - E;
    end
end
