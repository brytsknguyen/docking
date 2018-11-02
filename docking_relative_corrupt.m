close all;
clear all;

%% User-defined parameters

% Sampling period
tau = 0.1;

% Cycle
N = 96;

% Target's oscilation
M = 4096;

% Initial position in global frame
p0 = [-30; -30; 1.0];

% Target position in global
pT0 = [-5; 2; 10];

% Initial phase for the target
varphi = 0;

taralpha = 0;
tarbeta = 0;

% argument for f(d)
lambda = 100000000000000;

% Initial adaptive gain
gamma = 10;

% Proportional control gain
B = 5;

% Excitiation gain
alpha = 10;

% Maximum velocity
vM = 0.75;

% Stopping distance
d_stop = 0.0;

% Exciting magnitude
r1 = 0.5;
r2 = sqrt(3)/2;

% Target's periodic motion magnitude
tarPeriMag = diag([5, 5, 1.25]);

% Target's linear motion
vqL = 0.05*[-1; 1; 0.5];

% Initial phase of the autonomous PE process
psi = 0;

% Noises
% dmu = 0;
% dstd = 0;
% phimu = 0.0;
% phistd = 0.0;

dmu = 0;
dstd = 0.1;
phimu = 0.0;
phistd = 0.01;

% Relative docking position
qS = [-20, 5, 2]';

% Relative docking mode
redoc = 1;

%% Simulation
t = 0:tau:550;
% t = 0:tau:200;

K = size(t, 2);

G = gamma*eye(3, 3);

% Lyapunov function
V = zeros(1, K);

% Position in world frame
p = zeros(3, K);
p(:, 1) = p0;

% Relative position
q = zeros(3, K);

% Relative position estimate
qb = zeros(3, K);
qh = zeros(3, K);

qh(:, 1) = [0, 0, 0]';
qb(:, 1) = [0, 0, 0]';



% Velocity input
sigma = zeros(3, K);
rho = zeros(2, K);

u = zeros(3, K);
ub = zeros(3, K);

d = zeros(1, K);
dnoisy = zeros(1, K);
eps = zeros(1, K);
zeta = zeros(1, K);
phi = zeros(3, K);

% Dynamics of the autonomous process
D = [cos(2*pi/N), -sin(2*pi/N);...
     sin(2*pi/N), cos(2*pi/N)];


% Target's motion
pT = zeros(3, K);
pT(:, 1) = pT0 + tarPeriMag*[cos(taralpha); sin(taralpha); cos(tarbeta)];

figxyz = figure('position', [0, 0, 1000, 500], 'color', [1 1 1]);
hold on;
view([15, 9]);
grid on;
set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
ylabel('$y\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('$x\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
zlabel('$z\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
% daspect([1 1 1]);

uav_pos_gca = plot3(p(1, 1), p(2, 1), p(3, 1),...
                    '-.', 'color', [1, 0, 0], 'linewidth', 2);
                
uav_pos_est_gca = plot3(qh(1, 1) + pT(1, 1), qh(2, 1)+ pT(2, 1), qh(3, 1)+ pT(3, 1),...
                        '--', 'color', [0, 0, 1], 'linewidth', 2);

% tar_mark_gca = plot3(q(1, 1), q(2, 1), q(3, 1), 'o', 'markersize', 5, 'markerfacecolor', [1 0 0]);                   
if (norm(vqL) == 0) && (norm(tarPeriMag) == 0)
    tar_traj_gca = plot3(pT(1, 1), pT(2, 1), pT(3, 1), 'o', 'markersize', 10, 'markerfacecolor', [0 1 0]);
else
    tar_traj_gca = plot3(pT(1, 1), pT(2, 1), pT(3, 1), 'g', 'linewidth', 2);
end

% Mark the original locations of the trajectories
plot3(qh(1, 1) + pT(1, 1), qh(2, 1)+ pT(2, 1), qh(3, 1)+ pT(3, 1),...
                        'o', 'markerfacecolor', [0, 1, 0], 'markersize', 10);
plot3(p(1, 1), p(2, 1), p(3, 1),...
                    'o', 'markerfacecolor', [1, 0, 0], 'markersize', 10);
            
% reloc_gca = plot3(p(1, 1), p(2, 1), p(3, 1), '-',...
%                   'color', [1, 0, 0], 'linewidth', 2);
% reloc_est_gca = plot3(ph(1), ph(2), ph(3), '-',...
%                       'color', [0, 0, 1], 'linewidth', 2);

lhdxyz = legend('$p_k$', '$\hat{q}_k + p^*_k$', '$p^*_k$',...
                '$p^*_0$', '$p_0$');
            
set(lhdxyz, 'interpreter', 'latex', 'fontsize', 16,...
    'position', [0.1473 0.6791 0.1286 0.2862]);
% xlim([-40 0]);
% ylim([-30 11]);
daspect([1 1 1.5]);                     
% xlim([-5*Rc, 5*Rc]);
% ylim([-5*Rc, 5*Rc]);

k = 0;

reset = zeros(1, K);

premature_break = false;

percent = 0;

for k=1:K-1
    
    vqS = tarPeriMag*[cos(2*pi/M*(k-1) + taralpha);...
                      sin(2*pi/M*(k-1) + taralpha);...
                      cos(6*pi/M*(k-1) + tarbeta)];
	
	% Update target position
    pT(:, k+1) = pT0 + vqS + k*tau*vqL;
    
    % UWB distance measurement
    d(k) = norm(p(:, k) - pT(:, k));
    dnoisy(k) = d(k) + normrnd(dmu, dstd);
       
    if ceil(k/K*100) > percent
        percent = ceil(k/K*100);
        fprintf('Percent %d%%\n', percent);
    end
    
    % stop when distance to the target is below stopping distance
    if norm(p(:, k) - pT(:, k)) < d_stop
        fprintf('Stopping time: %f\n', k*tau);
        premature_break = true;
        break;
    end
    
    % Generate velocity setpoint
    if k == 1
        rho(:, k) = [cos(varphi); sin(varphi)];
%         sigma(:, k) = [r1*(4*rho(1, k)^3 - 3*rho(1, k));...
%                        r1*(3*rho(2, k) - 4*rho(2, k)^3);...
%                        r2*rho(1, k)];

        sigma(:, k) = [r1*rho(1, k);...
                       r1*rho(2, k);...
                       r2*(4*rho(1, k)^3 - 3*rho(1, k))];

        
        q(:, 1) = p0(:, 1) - pT(:, 1);
        
%         fd = 1/lambda*(1 - exp(-lambda*dnoisy));
        fd = lambda/max([dnoisy(k) lambda])*dnoisy(k);
        
        if redoc
            u(:, k) = - B*(qh(:, k) - qS) + alpha*sigma(:, k);
        else
            u(:, k) = - B*qh(:, k) + fd*alpha*sigma(:, k);
        end
        
    else
        
        rho(:, k) = D*rho(:, k-1);
        sigma(:, k) = [r1*(4*rho(1, k)^3 - 3*rho(1, k));...
                       r1*(3*rho(2, k) - 4*rho(2, k)^3);...
                       r2*rho(1, k)];
        
        phi(:, k) = [p(:, k) - p(:, k-1)] + [normrnd(phimu, phistd);...
                                             normrnd(phimu, phistd);...
                                             normrnd(phimu, phistd)];
        
        zeta(k) = 0.5*(dnoisy(k)^2 - dnoisy(k-1)^2 - phi(:, k)'*phi(:, k));

        eps(k) = zeta(k) - phi(:, k)'*qh(:, k-1);
            
        qb(:, k) = qh(:, k-1) + phi(:, k) + G*phi(:, k)*eps(k);
        qh(:, k) = dnoisy(k)/max([norm(qb(:, k)), dnoisy(k)])*qb(:, k);
%         qh(:, k) = qb;
        if norm(qh(:, k) - qb(:, k)) ~= 0
            fprintf('projection:%d\n', norm(qh(:, k) - qb(:, k)));
%             while(true)
%             end
        end

%         fd = 1/lambda*(1 - exp(-lambda*dnoisy(k)));

        fd = lambda/max([dnoisy(k) lambda])*dnoisy(k);
        
        if redoc
            u(:, k) = - B*(qh(:, k) - qS) + alpha*sigma(:, k);
        else
            u(:, k) = - B*qh(:, k) + fd*alpha*sigma(:, k);
        end
        
    end
    
    if norm(u(:, k)) > vM 
        ub(:, k) = u(:, k)/norm(u(:, k))*vM;
    else
        ub(:, k) = u(:, k);
    end
        
    % Update the absolute position
    p(:, k+1) = p(:, k) + tau*ub(:, k);
    
    % Update the relative position
    q(:, k+1) = p(:, k+1) - pT(:, k+1);
    
    % Update the Lyapunov function
    qt = qh(:, k) - q(:, k);
    V(k) = qt'*G^(-1)*qt;
    
    % Illustrate the latest update
    
%     % Target location and trajectory
%     set(tar_mark_gca, 'xdata', pT(1, k+1),...
%                       'ydata', pT(2, k+1),...
%                       'zdata', pT(3, k+1));
%                   
%     set(tar_traj_gca, 'xdata', pT(1, 1:k+1),...
%                       'ydata', pT(2, 1:k+1),...
%                       'zdata', pT(3, 1:k+1));
%                   
%     % Target location estimate and UAV location
%     set(tar_est_traj_gca, 'xdata', pTh(1, 1:k),...
%                           'ydata', pTh(2, 1:k),...
%                           'zdata', pTh(3, 1:k));
%                       
%     set(uav_pos_traj_gca, 'xdata', p(1, 1:k+1),...
%                           'ydata', p(2, 1:k+1),...
%                           'zdata', p(3, 1:k+1));
%     pause(min([0.01, tau]));
end

t(k:end) = [];
ub(:, k:end) = [];
u(:, k:end) = [];
p(:, k:end) = [];
pT(:, k:end) = [];

q(:, k:end) = [];

qb(:, k:end) = [];
qh(:, k:end) = [];
phi(:, k:end) = [];
zeta(:, k:end) = [];
eps(:, k:end) = [];
d(k:end) = [];
dnoisy(k:end) = [];
V(k:end) = [];
reset(k:end) = [];
sigma(:, k:end) = [];

% set(tar_mark_gca, 'xdata', q(1, k-1),...
%     'ydata', q(2, k-1),...
%     'zdata', q(3, k-1));

set(tar_traj_gca, 'xdata', pT(1, :),...
                  'ydata', pT(2, :),...
                  'zdata', pT(3, :));

% UAV absolute position
set(uav_pos_gca, 'xdata', p(1, :),...
                 'ydata', p(2, :),...
                 'zdata', p(3, :));
              
% UAV estimated position
set(uav_pos_est_gca, 'xdata', qh(1, :) + pT(1, :),...
                     'ydata', qh(2, :) + pT(2, :),...
                     'zdata', qh(3, :) + pT(3, :));

% % Relative position and its estimate
% set(reloc_est_gca, 'xdata', ph(1, :),...
%                    'ydata', ph(2, :),...
%                    'zdata', ph(3, :));
% 
% set(reloc_gca, 'xdata', p(1, :),...
%                'ydata', p(2, :),...
%                'zdata', p(3, :));


% figure('position', [1920, 1080/4, 1080, 480], 'color', [1 1 1]);
% hold on;
% plot(t, ub(1, :));
% plot(t, ub(2, :));
% title('velocity');

%% xyz plot
if ~redoc
    figure('position', [100, 0, 1080, 720], 'color', [1 1 1]);
    subplot(3, 1, 1);
    hold on;
    plot(t, q(1, :), '-.r', 'linewidth', 2);
    plot(t, qh(1, :), '--b', 'linewidth', 2);
    plot(t, qh(1, :) - q(1, :), 'g', 'linewidth', 2);
    % plot(t, pw(1, :), 'b', 'linewidth', 2);
    ylabel('$x\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
    % xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
    lhdx = legend('$q_k$', '$\hat{q}_k$', '$\tilde{q}_k$');
    set(lhdx, 'interpreter', 'latex', 'fontsize', 20, 'position', [0.8737, 0.7, 0.0803, 0.1162]);
    xlim([0, t(end)]);
    % ylim([-35, 10]);
    grid on;

    subplot(3, 1, 2);
    grid on;
    hold on;
    plot(t, q(2, :), '-.r', 'linewidth', 2);
    plot(t, qh(2, :), '--b', 'linewidth', 2);
    plot(t, qh(2, :) - q(2, :), 'g', 'linewidth', 2);
    ylabel('$y\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
    % xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
    % lhdy = legend('$p_{T, y}$', '$\hat{q}_{y}$', '$p_{y}$');
    % set(lhdy, 'interpreter', 'latex', 'fontsize', 20);
    xlim([0, t(end)]);
    % ylim([-35, 10]);
    grid on;

    subplot(3, 1, 3);
    grid on;
    hold on;
    plot(t, q(3, :), '-.r', 'linewidth', 2);
    plot(t, qh(3, :), '--b', 'linewidth', 2);
    plot(t, qh(3, :) - q(3, :), 'g', 'linewidth', 2);
    ylabel('$z\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
    xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
    % lhdz = legend('$p_{T, z}$', '$\hat{q}_{z}$', '$p_{z}$');
    % set(lhdz, 'interpreter', 'latex', 'fontsize', 20);
    xlim([0, t(end)]);
    % ylim([0, 12]);
    grid on;
    tightfig(gcf);
else
    %% redoc xyz
    figure('position', [100, 0, 1080, 720], 'color', [1 1 1]);
    subplot(3, 1, 1);
    hold on;
    plot(t, q(1, :) - qS(1), '-.r', 'linewidth', 2);
    plot(t, qh(1, :) - qS(1), '--b', 'linewidth', 2);
    plot(t, qh(1, :) - q(1, :), 'g', 'linewidth', 2);
    % plot(t, pw(1, :), 'b', 'linewidth', 2);
    ylabel('$x\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
    % xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
    xlim([0, t(end)]);
    % ylim([-35, 10]);
    grid on;

    subplot(3, 1, 2);
    grid on;
    hold on;
    plot(t, q(2, :) - qS(2), '-.r', 'linewidth', 2);
    plot(t, qh(2, :) - qS(2), '--b', 'linewidth', 2);
    plot(t, qh(2, :) - q(2, :), 'g', 'linewidth', 2);
    ylabel('$y\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
    % xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
    % lhdy = legend('$p_{T, y}$', '$\hat{q}_{y}$', '$p_{y}$');
    % set(lhdy, 'interpreter', 'latex', 'fontsize', 20);
    xlim([0, t(end)]);
    % ylim([-35, 10]);
    lhdx = legend('$\bar{q}_k$', '$\hat{q}_k - q^*$', '$\tilde{q}_k$');
    set(lhdx, 'interpreter', 'latex', 'fontsize', 20, 'position', [0.8392 0.6578 0.1407 0.1500]);
    grid on;

    subplot(3, 1, 3);
    grid on;
    hold on;
    plot(t, q(3, :) - qS(3), '-.r', 'linewidth', 2);
    plot(t, qh(3, :) - qS(3), '--b', 'linewidth', 2);
    plot(t, qh(3, :) - q(3, :), 'g', 'linewidth', 2);
    ylabel('$z\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
    xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
    % lhdz = legend('$p_{T, z}$', '$\hat{q}_{z}$', '$p_{z}$');
    % set(lhdz, 'interpreter', 'latex', 'fontsize', 20);
    xlim([0, t(end)]);
    % ylim([0, 12]);
    grid on;
    tightfig(gcf);
end
% uav_pos_traj_gca = plot3(p(1, 1), p(2, 1), p(3, 1), '-',...
%                          'color', [0, 0.5, 0]);
% uav_pos_est_traj_gca = plot3(p(1), p(2), p(3), '-',...
%                             'color', [0.5, 0.5, 0]);

%% redoc error
qbar = [q(1, :) - qS(1);...
        q(2, :) - qS(2);...
        q(3, :) - qS(3)];

qbar_nrm = sqrt(dot(qbar, qbar, 1));

alphabarS = alpha/B*(1 + abs(1 - tau*B))/(2-tau*B);

figure('position', [200, 0, 800, 300], 'color', [1 1 1]);
hold on;
plot(t, qbar_nrm, 'b', 'linewidth', 2);
plot([t(1),t(end)], [alphabarS, alphabarS], 'r', 'linewidth', 2);
% title('Distance');
grid on;
set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
% ylabel('$||\bar{q}_k||\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$[m]$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
xlim([0, t(end)]);
lhdredocnorm = legend('$\|\bar{q}_k\|$', '$\bar{\alpha}$');
set(lhdredocnorm, 'interpreter', 'latex',...
                  'fontsize', 20, 'Position', [0.3786 0.6846 0.1488 0.2350]);

axes('Position', [0.5724 0.5065 0.3318 0.4231])
box on
hold on;
plot(t(round(5*end/6):end), qbar_nrm(round(5*end/6):end), 'b', 'linewidth', 2);
plot([t(round(5*end/6)),t(end)], [alphabarS, alphabarS], 'r', 'linewidth', 2);
ylim([0, 8]);
xlim([t(round(5*end/6)), t(end)]);
set(gca, 'Fontname', 'cambria math', 'fontsize', 16);
grid on;
tightfig(gcf);


%% Plot distance
% find the last point d is outside 1m
KK = size(t, 2);
ek = 1;
for k =1:KK
    if d(k) > 3
        ek = k;
    end
end
fprintf('entry time: %f, %f\n', t(ek), d(ek));

figure('position', [200, 0, 800, 300], 'color', [1 1 1]);
hold on;
plot(t, d, 'linewidth', 2);
% title('Distance');
grid on;
set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
ylabel('$d_k\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
xlim([0, t(end)]);
tightfig(gcf);


figure('position', [300, 0, 1080, 480], 'color', [1 1 1]);
hold on;
dV = V(2:end) - V(1:end-1);
plot(t, V, 'linewidth', 2);
plot(t(2:end-1), dV(2:end), 'linewidth', 2);
% title('V and dV');
set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
xlabel('$\mathrm{Time}\ [s]$', 'interpreter', 'latex', 'fontsize', 20);
lhdV = legend('$V_k$', '$\Delta V_k$');
set(lhdV, 'interpreter', 'latex', 'fontsize', 20,...
    'position', [0.7771, 0.6117, 0.0786, 0.1613]);
grid on;
tightfig(gcf);


% figure;
% hold on;
% plot3(sigma(1, :), sigma(2, :), sigma(3, :));

% figure;
% hold on;
% plot(t, reset*10);
% title('reset times');

% Some pausing to make sure the figure's content is settled before we
% tighten it
pause(0.01);
tightfig(figxyz);

sigcol = hsv(N);
figure('position', [200, 0, 600, 900], 'color', [1 1 1]);
hold on;
% scatter3(sigma(1, 1:N), sigma(2, 1:N), sigma(3, 1:N), 20, [1:N]/N, 'filled');
plot3(sigma(1, :), sigma(2, :), sigma(3, :),...
      'linewidth', 0.5, 'color', 'r');
for k = 1:N
  quiver3(0, 0, 0, 1.1*sigma(1, k), 1.1*sigma(2, k), 1.1*sigma(3, k),...
      'linewidth', 1.5, 'color', sigcol(k, :));
end
set(gca, 'Fontname', 'cambria math', 'fontsize', 20);
xlabel('$x\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
zlabel('$z\ [m]$', 'interpreter', 'latex', 'fontsize', 20);
grid on;
daspect([1 1 1]);
view([69 37]);
pause(0.01);
tightfig(gcf);