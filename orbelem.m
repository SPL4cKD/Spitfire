function [h, e, i, w, RAAN, theta] = orbelem(r, v, mu)
%Orbit Elements, Finds a sattellite's orbit elements from a state vector
%

h = cross(r,v)';
n = cross([0 0 1],h)';
e = cross(v,h)'/mu - r/norm(r);
i = acosd(h(3)/norm(h));

if e(3) < 0
    w = rad2deg(2*pi-acos(dot(n,e)/(norm(n)*norm(e))));
else
    w = acosd(dot(n,e)/(norm(n)*norm(e)));
end

if n(2) < 0
    RAAN = rad2deg(2*pi -acos(n(1)/norm(n)));
else
    RAAN = acosd(n(1)/norm(n));
end

if dot(r,v) < 0
    theta = rad2deg(2*pi - acos(dot(e,r)/(norm(e)*norm(r))));
else
    theta = acosd(dot(e,r)/(norm(e)*norm(r)));
end

end

