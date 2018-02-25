function [ ] = orbelemP(h, e, i, w, RAAN, theta )
%Orb Element Print, Prints the orbit elements in a nice format
%

fprintf('Angular Momentum = %g  [km^2/s]\n',norm(h));
fprintf('Eccentricity     = %g\n',norm(e));
fprintf('Inclination      = %g  [deg]\n',i);
fprintf('Angle of Perigee = %g  [deg]\n',w);
fprintf('RAAN             = %g  [deg]\n',RAAN);
fprintf('True Anomaly     = %g  [deg]\n\n',theta);

end

