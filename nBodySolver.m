function [t,x,dx,ddx] = nBodySolver(p0, v0, mu, tt)

leng = length(mu);
x123 = repmat(1:leng,1,leng);
x1214 = reshape(1:leng*(leng-1),leng-1,leng)';
x23 = x123(logical(ones(leng) - eye(leng)));
x11 = reshape(repmat(1:leng,leng-1,1),1,leng*(leng-1));

% [t,dz]   = ode113(@nBodyFunc,tt,[p0;v0]);
% [t,x,dx]   = rkn86(@nBodyFunc,tt(1),tt(end),p0,v0);
% [t,dz]   = ode23(@nBodyFunc,tt,[p0;v0]);
[t,x,dx]   = rkn1210(@nBodyFunc2,tt,p0,v0);


% x  = dz(:,1:3*leng);
% dx = dz(:,3*leng+1:end);

xSize = size(x);
ddx = zeros(xSize);
for i = 1:xSize(1)
    xv = [x(i,:)';dx(i,:)'];
    pos = reshape(xv(1:3*leng),3,leng);
    r = pos(:,x23) - pos(:,x11);
    r3 = diag(1./sqrt(sum(r.^2)).^3);
    gm = diag(mu(x23));
    dpos = r*r3*gm;
    ddx(i,:) = sum(reshape(dpos(:,x1214(:)),3*leng,leng-1),2);
end

    function [ddpos] = nBodyFunc2(t,xv)
        pos = reshape(xv(1:3*leng),3,leng);
        r = pos(:,x23) - pos(:,x11);
        r3 = diag(1./sqrt(sum(r.^2)).^3);
        gm = diag(mu(x23));
        dpos = r*r3*gm;
        ddpos = sum(reshape(dpos(:,x1214(:)),3*leng,leng-1),2);
    end

    function [dz] = nBodyFunc(t,xv)
        pos = reshape(xv(1:3*leng),3,leng);
        r = pos(:,x23) - pos(:,x11);
        r3 = diag(1./sqrt(sum(r.^2)).^3);
        gm = diag(mu(x23));
        dpos = r*r3*gm;
        ddpos = sum(reshape(dpos(:,x1214(:)),3*leng,leng-1),2);
        dz = [xv(3*leng+1:end);ddpos];
    end
end

