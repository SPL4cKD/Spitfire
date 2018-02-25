%%
%

% while 1
clear all;
close all;
clc;

%% Input
%

% iterations and speed
iter   = 1e4;
spdn   = 2;

% Time and Iterations
yrs    = 0;        % years-
days   = 10;        % days
hrs    = 0;        % hours

ApNum = 3;

% while 1
%     try
%         ApNum = input('Which Simulation: '); % rePlay input
%     catch ME
%         fprintf(2,'%s\n', ME.message);
%         continue
%     end
%     if isempty(ApNum) % Usable value test
%         fprintf(2,'Error: input not recognized.\n'); % Error notification
%         continue
%     elseif ApNum >= 10 && ApNum <= 12 && rem(ApNum,1) == 0 || ApNum == 5 || ApNum == 6 || ApNum == 0 || ApNum == 1 || ApNum == 3
%         break
%     else
%         fprintf(2,'Error: input not recognized.\n'); % Error notification
%         continue
%     end
% end

ref_Frame = 1;

% if ApNum ~= 0
%     while 1
%         try
%             ref_Frame = input('Frame [Inertial(1), Rotational(2)]: '); % rePlay input
%         catch ME
%             fprintf(2,'%s\n', ME.message);
%             continue
%         end
%         if isempty(ref_Frame) % Usable value test
%             fprintf(2,'Error: input not recognized.\n'); % Error notificati1on
%             continue
%         elseif ref_Frame == 1 || ref_Frame == 2
%             break
%         else
%             fprintf(2,'Error: input not recognized.\n'); % Error notification
%             continue
%         end
%     end
% else
%     ref_Frame = 1;
% end



%% Initial State Vectors and Masses
%

if ApNum >= 10 && ApNum <= 12 && rem(ApNum,1) == 0;
    [p0, v0, mu, n, scale, cA, bt] = ApolloCoords(ApNum);
    
%     if ApNum == 10
%         [p0, v0, mu, n, scale, cA ] = bodyAdd(p0,v0,mu,n,scale,cA);
%     end
    pFrame = 4:6;
    [ap10x] = Ap10(ApNum);
elseif ApNum == 5;
    [p0, v0, mu, n, scale, cA, bt] = LagrangeP(ApNum);
    pFrame = 4:6;
elseif ApNum == 6;
    [p0, v0, mu, n, scale, cA, bt] = CassiniCoord;
% elseif ApNum == 1;
%     [p0, v0, mu, n, scale, cA, bt, pMinus, M] = ringBS;
elseif ApNum == 0;
    [p0, v0, mu, n, scale, cA, bt] = StateVecInit;
    pFrame = 10:12;
elseif ApNum == 3;
    [p0, v0, mu, n, scale, cA, bt, graAr] = FallSolEr;
    [p0e, v0e, mue, ne, scalee, cAe, bte, graAre] = FallNoP;
    pFrame = 1:3;
end

%% Time Parameters
%


% Time span array calculation
time    = [0, ((yrs*365.25+days)*24+hrs)*3600]; % Convert to seconds                    % Time of 1st Burn
% time = ((yrs*365.25+days)*24+hrs)*3600;
[tt, tn, solv2] = SatBurn(time,bt,iter);
% tt = linspace(0,time,iter);

% if time <= b1Time || ApNum ~= 11
%     solv2 = 1;
%     tt1       = linspace(0, time, iter);
% else
%     time1 = b1Time;
%     time2 = time;
%     iter1 = round(b1Time/time*iter);
%     iter2 = round((1-b1Time/time)*iter);
%     solv2 = 2;
%     tt1   = linspace(0, time1, iter1);
%     tt2   = linspace(time1, time2, iter2);
% end

spd       = round(spdn*iter/1e3);
satDotDiv = time(end)/4.32e7;
dotNS     = [30 100];               % [num of dots, spacing of dots]

%% N-Body Function for Numerical Integration Solver
%
tic
fprintf('\n  working ...\n\n');


% Numerical Solver Function
[t1, x1, dx1, ddx1] = nBodySolver( p0, v0, mu, tt{1} );
[t1e, x1e, dx1e, ddx1e] = nBodySolver( p0e, v0e, mue, tt{1} );

t = t1;
x = x1;
dx = dx1;
ddx = ddx1;

te = t1e;
xe = x1e;
dxe = dx1e;
ddxe = ddx1e;
solv2 = 0;
if solv2 == 2
    p02 = x1(end,:)';
    v02 = dx1(end,:)';
    p02e = x1e(end,:)';
    v02e = dx1e(end,:)';
    if ApNum == 11
        v02(7:9) = v02(7:9)*0.2139985809;
    end
    
    [t2, x2, dx2, ddx2] = nBodySolver( p02, v02, mu, tt{2});
    [t2e, x2e, dx2e, ddx2e] = nBodySolver( p02e, v02e, mue, tt{2});
    
    t = vertcat(t1(1:end-1,:),t2);
    x = vertcat(x1(1:end-1,:),x2);
    dx = vertcat(dx1(1:end-1,:),dx2);
    ddx = vertcat(ddx1(1:end-1,:),ddx2);
    
    te = vertcat(t1e(1:end-1,:),t2e);
    xe = vertcat(x1e(1:end-1,:),x2e);
    dxe = vertcat(dx1e(1:end-1,:),dx2e);
    ddxe = vertcat(ddx1e(1:end-1,:),ddx2e);
end



figN = 1;
if ref_Frame == 2
    
    
    xR = zeros(size(x));
    dxR = zeros(size(x));
    ddxR = zeros(size(x));
    Arot = vrrotvec2mat(vrrotvec(x(1,pFrame),[-1 0 0]));
    for i = 1:size(x,2)/3
        
        ip = 1+(i-1)*3;
        fp = 3+(i-1)*3;
        for j = 1:spd:size(x,1)
            
            v1x = x(j,pFrame);
            vec = x(1,pFrame);
            
            Arot2 = Arot*vrrotvec2mat(vrrotvec(v1x,vec));
            xR(j,ip:fp) = Arot2*x(j,ip:fp)';
            dxR(j,ip:fp) = Arot2*dx(j,ip:fp)';
            ddxR(j,ip:fp) = Arot2*ddx(j,ip:fp)';
            
        end
    end
    
    Arot = vrrotvec2mat(vrrotvec([0, xR(1,8:9)],[0 1 0]));
    for i = 1:size(xR,2)/3
        
        ip = 1+(i-1)*3;
        fp = 3+(i-1)*3;
        for j = 1:spd:size(xR,1)
            
            v1x = xR(j,pFrame);
            vec = xR(1,pFrame);
            
            Arot2 = Arot*vrrotvec2mat(vrrotvec(v1x,vec));
            xR(j,ip:fp) = Arot2*xR(j,ip:fp)';
            dxR(j,ip:fp) = Arot2*dxR(j,ip:fp)';
            ddxR(j,ip:fp) = Arot2*ddxR(j,ip:fp)';
            
        end
    end
    
    x = xR;
    dx = dxR;
    ddx = ddxR;
    
    subplot(3,1,1)
    plot(x(1:spd:end,7))
    subplot(3,1,2)
    plot(x(1:spd:end,8))
    subplot(3,1,3)
    plot(x(1:spd:end,9))
    
    set(gcf,'position',[649 70 1257 918]);
    figN = 1;
end

toc
%% Graphical Setup

% Change of reference frame
efre = 1:3;
efr = 10:12;
mfr = 4:6;
sfr = 7:9;

% satfr = 1+3*12:3+3*11;
mDiff = x(:,efr) - repmat(x(1,efr),size(x,1),1);

% satfre = 1+3*11:3+3*11;
mDiffe = xe(:,efre) - repmat(xe(1,efre),size(xe,1),1);

x = x - repmat(mDiff,1,size(x,2)/3);
xe = xe - repmat(mDiffe,1,size(xe,2)/3);


fprintf('\n___Initial Orbit Parameteres of ISS___\n')
fprintf('ISS (initial) with Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(x(1,1+11*3:3+11*3)', dx(1,1+11*3:3+11*3)', mu(4));
orbelemP(h, e, i, w, RAAN, theta )
fprintf('ISS (initial) w/o Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(xe(1,4:6)', dxe(1,4:6)', mue(1));
orbelemP(h, e, i, w, RAAN, theta )
fprintf('ISS (final) with Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(x(end,1+11*3:3+11*3)', dx(end,1+11*3:3+11*3)', mu(4));
orbelemP(h, e, i, w, RAAN, theta )
fprintf('ISS (final) w/o Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(xe(end,4:6)', dxe(end,4:6)', mue(1));
orbelemP(h, e, i, w, RAAN, theta )
fprintf('\n\n')

fprintf('___Initial Orbit Parameteres of Chandra___\n')
fprintf('Chandra (initial) with Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(x(1,1+12*3:3+12*3)', dx(1,1+12*3:3+12*3)', mu(4));
orbelemP(h, e, i, w, RAAN, theta )
fprintf('Chandra (initial) w/o Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(xe(1,7:9)', dxe(1,7:9)', mue(1));
orbelemP(h, e, i, w, RAAN, theta )

fprintf('Chandra (final) with Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(x(end,1+12*3:3+12*3)', dx(end,1+12*3:3+12*3)', mu(4));
orbelemP(h, e, i, w, RAAN, theta )
fprintf('Chandra (final) w/o Planet Perturbations\n')
[h, e, i, w, RAAN, theta] = orbelem(xe(end,7:9)', dxe(end,7:9)', mue(1));
orbelemP(h, e, i, w, RAAN, theta )




ylimMin = round(min(min(x(:,34:36))));
ylimMax = round(max(max(x(:,34:36))));
ylimMin = round(ylimMin,-(numel(num2str(abs(ylimMin)))-1));
ylimMax = round(ylimMax,-(numel(num2str(abs(ylimMax)))-1));
ylimRange = [ylimMin, ylimMax];

set(0,'defaultfigurecolor', [1 1 1])
figure('units','normalized','outerposition',[0 0 1 1])
set(gca,'Color', [1 1 1])

subplot(3,1,1)
plot(t, x(:,1+11*3)); hold on
plot(te, xe(:,4),'--')
ylabel('X [km]')
xlabel('Time [sec]')
ylim(ylimRange*1.5)
title('ISS XYZ Position vs Time')
legend('Perturb','w/o perturb','location','best')

subplot(3,1,2)
plot(t, x(:,2+11*3)); hold on
plot(te, xe(:,5),'--')
ylabel('Y [km]')
xlabel('Time [sec]')
ylim(ylimRange*1.5)
legend('Perturb','w/o perturb','location','best')

subplot(3,1,3)
plot(t, x(:,3+11*3)); hold on
plot(te, xe(:,6),'--')
ylabel('Z [km]')
xlabel('Time [sec]')
ylim(ylimRange*1.5)
legend('Perturb','w/o perturb','location','best')








ylimMin = round(min(min(x(:,37:39))));
ylimMax = round(max(max(x(:,37:39))));
ylimMin = round(ylimMin,-(numel(num2str(abs(ylimMin)))-1));
ylimMax = round(ylimMax,-(numel(num2str(abs(ylimMax)))-1));
ylimRange = [ylimMin, ylimMax];

set(0,'defaultfigurecolor', [1 1 1])
figure('units','normalized','outerposition',[0 0 1 1])
set(gca,'Color', [1 1 1])

subplot(3,1,1)
plot(t, x(:,1+12*3)); hold on
plot(te, xe(:,7),'--')
ylabel('X [km]')
xlabel('Time [sec]')
ylim(ylimRange*1.3)
title('Chandra XYZ Position vs Time')
legend('Perturb','w/o perturb','location','best')

subplot(3,1,2)
plot(t, x(:,2+12*3)); hold on
plot(te, xe(:,8),'--')
ylabel('Y [km]')
xlabel('Time [sec]')
ylim(ylimRange*1.3)
legend('Perturb','w/o perturb','location','best')

subplot(3,1,3)
plot(t, x(:,3+12*3)); hold on
plot(te, xe(:,9),'--')
ylabel('Z [km]')
xlabel('Time [sec]')
ylim(ylimRange*1.3)
legend('Perturb','w/o perturb','location','best')
figN = 1;

while 1
    while 1
        try
            anim = input('Animation [Yes(1), No(2)]: '); % rePlay input
        catch ME
            fprintf(2,'%s\n', ME.message);
            continue
        end
        if isempty(anim) % Usable value test
            fprintf(2,'Error: input not recognized.\n'); % Error notificati1on
            continue
        elseif anim == 1 || anim == 2
            break
        else
            fprintf(2,'Error: input not recognized.\n'); % Error notification
            continue
        end
    end
    if anim == 1
        plst = 1;
    elseif anim == 2
        plst = size(x,1);
    end
    
    wmax = max(max(x(:,:))) * 1.1;
    wmin = min(min(x(:,:))) * 2;
    ax_range = [wmin,-wmin, wmin,-wmin, wmin,-wmin];
    
    if figN == 1
        
        pRadius = 6378.1; % [km] central planet radius
        
        [u1,n1,w1] = sphere;
        us = u1;
        vs = n1;
        ws = w1;
        
        set(0,'defaultfigurecolor', [0 0 0])
        figure('units','normalized','outerposition',[0 0 1 1])
%         set(gcf,'position',[649 70 1257 918]);
%         set(gcf,'position',[-2559 167 2560 1006]);
        axis equal
        %             axis(ax_range,'square')
        grid on
        set(gca,'Color', [0 0 0])
        ax2 = gca;
        alpha(1)
        
        xlabel('X [km]');
        ylabel('Y [km]');
        zlabel('Z [km]');
        rotate3d on
        figN = 0;
    end
    %% Animation Loop for Orbit Visualization
    %
    
    
    pColor = SetColor(n);
    plotStep = plst:spd:size(x,1);
    plotStepe = plst:spd:size(xe,1);
    st2 = round(dotNS(2)/satDotDiv);
    k2 = repmat(1:st2,1,ceil(dotNS(2)*dotNS(1)/st2))';
    h1 = zeros(plotStep(end),n);
    h2 = zeros(plotStep(end),n);
    h1e = zeros(plotStepe(end),ne);
    h2e = zeros(plotStepe(end),ne);
    
    if ApNum >= 10 && ApNum <= 12 && ref_Frame == 1
        h3 = plot3(ap10x(:,1), ap10x(:,2), ap10x(:,3), '-', 'color', pColor(4,1:3)); hold on
    end
    
    for i = plotStep
        
        [k] = TrailSpacing(i,k2,dotNS);
        
        if n > 1
            for j = graAr
                
                h1(i,j) = surf(us*scale(j) + x(i,1+(j-1)*3), vs*scale(j) + x(i,2+(j-1)*3), ws*scale(j) + x(i,3+(j-1)*3), 'FaceColor', 'none', 'EdgeColor', pColor(cA(j),1:3)); hold on
                % h2(i,j) = plot3(x(k:st2:i,1+(j-1)*3), x(k:st2:i,2+(j-1)*3), x(k:st2:i,3+(j-1)*3), '.', 'markeredgecolor', pColor(cA(j),1:3), 'markers', .01); hold on
                h2(i,j) = plot3(x(1:spd:i,1+(j-1)*3), x(1:spd:i,2+(j-1)*3), x(1:spd:i,3+(j-1)*3), '.', 'markeredgecolor', pColor(cA(j),1:3), 'markers', .01); hold on
            end
        end
        
        for j = graAre
            
            h1e(i,j) = surf(us*scalee(j) + xe(i,1+(j-1)*3), vs*scalee(j) + xe(i,2+(j-1)*3), ws*scalee(j) + xe(i,3+(j-1)*3), 'FaceColor', 'none', 'EdgeColor', pColor(cAe(j),1:3)); hold on
            % h2(i,j) = plot3(x(k:st2:i,1+(j-1)*3), x(k:st2:i,2+(j-1)*3), x(k:st2:i,3+(j-1)*3), '.', 'markeredgecolor', pColor(cA(j),1:3), 'markers', .01); hold on
            h2e(i,j) = plot3(xe(1:spd:i,1+(j-1)*3), xe(1:spd:i,2+(j-1)*3), xe(1:spd:i,3+(j-1)*3), '.', 'markeredgecolor', pColor(cAe(j),1:3), 'markers', .01); hold on
        end
        % Gravitational Potential Surface
%                 [X,Y] = meshgrid(-wmin:wmin/100:wmin, -wmin:wmin/100:wmin);
%                 Z = -(mu(1)./(sqrt((X-x(i,1)).^2+(Y-x(i,2)).^2)).^1 + mu(2)./(sqrt((X-x(i,4)).^2+(Y-x(i,5)).^2).^1))*1e5-1e4;
%                 for zi = 1:size(Z,1)
%                     for zj = 1:size(Z,2)
%                         if Z(zi,zj) < -3e5
%                             Z(zi,zj) = -3e5;
%                         end
%                     end
%                 end
%                 h4 = surf(X,Y,Z,'facecolor',[.3 .3 .3],'FaceAlpha',0.1,'EdgeColor',[.1 .1 .1]);
        
        grid on
        axis equal
        %             axis(ax_range,'square')
        set(gca,'Color',[0 0 0])
        ax2.XColor = [.3, .3, .3];
        ax2.YColor = [.3, .3, .3];
        ax2.ZColor = [.3, .3, .3];
        ax2.GridAlpha = 0.15;
        ax2 = gca;
        ax2.GridColor = [1 1 1];
        alpha(1)
        xlabel('X [km]');
        ylabel('Y [km]');
        zlabel('Z [km]');
         
        pause(1e-20)
        
        if i ~= plotStep(end)
            delete(h1(i,:))
            delete(h2(i,:))
            delete(h1e(i,:))
            delete(h2e(i,:))
%             delete(h3)
%             delete(h4)
        end
    end
    
    
    % Requests user input to stop or run the animation again
    while 1
        try
            rePlay = input('\nWould you like to graph again? \n Yes = 1, No = 0: '); % rePlay input
        catch ME
            fprintf(2,'%s\n', ME.message);
            continue
        end
        if isempty(rePlay) % Usable value test
            fprintf(2,'Error: input not recognized.\n'); % Error notification
            continue
        elseif rePlay == 1
            fprintf('\n')
            delete(h1(i,:))
            delete(h2(i,:))
            delete(h1e(i,:))
            delete(h2e(i,:))
%             delete(h3)
%             delete(h4)
            break
        elseif rePlay == 0
            break
        else
            fprintf(2, 'Error: Value was not 1 or 0.\n');
            continue
        end
    end
    if rePlay == 0
        break
    end
%     set(0,'defaultfigurecolor', [1 1 1])
end
% end