function ds = cr3bp_ode(~, s, mu)
    x = s(1); y = s(2); z = s(3);
    vx = s(4); vy = s(5); vz = s(6);
    
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    
    ax = x + 2*vy - (1-mu)*(x + mu)/r1^3 - mu*(x - 1 + mu)/r2^3;
    ay = y - 2*vx - (1-mu)*y/r1^3 - mu*y/r2^3;
    az = -(1-mu)*z/r1^3 - mu*z/r2^3;
    
    ds = [vx; vy; vz; ax; ay; az];
end