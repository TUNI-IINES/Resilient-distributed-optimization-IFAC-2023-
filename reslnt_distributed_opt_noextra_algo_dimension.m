


%==========================================================
%|     Continuous-Time Distributed Optimization with attack
%==========================================================

clear all;
close all;
clc;
n = 8;
x_ini = [1.1 2 2.3 2 1.4 2.5 1 2 1.7 2.8 5.6 2 6 2 6 1]';

z_ini = [1 2 2 2 1.4 2.7 1 2.6 1.5 2.3 5.1 2.1 6.5 2.7 3 2]';
d_ini = rand(2*n,1);
chi_ini = [ x_ini ; z_ini ; d_ini];

tspan = [0 4];
options = odeset('RelTol',1e-7);
[t, chi] = ode45(@rslnt_iman, tspan, chi_ini, options);

opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 17;
opts.height     = 5;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% create new figure
fig = figure; clf

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;
axis tight
% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
% export to png
fig.PaperPositionMode   = 'auto';
plot(t, chi(:, 1:16),'linewidth',1.5);


xlabel('Time (Seconds)','FontSize',16);
legend('$\boldsymbol{x}_i$', '$\boldsymbol{z}_i$', '$\text{for all} \; i \in \mathcal V $','Interpreter','latex','FontSize',16)
%ylabel('$x_i $ for all $ i \in \mathcal V $' , 'Interpreter','latex','FontSize',16);
% save2pdf('privacy2_low_beta2',fig,600);
save2pdf('distributed_opt_highbeta',fig,600);
%%%%%%%%%%%%%%%%%%%%
function chi_dot = rslnt_iman(t, chi)

beta = 300;





syms x1 x2
fx1 = (x1 + x2)^2;
fx2 = (2*x1 + 3*x2 + 1)^2;
fx3 = (x1 + 2*x2 - 3)^2;
fx4 = (3*x1 + x2 + 5)^2;
fx5 = (4*x1 + x2 + 1)^2;
fx6 = (x1 - 5*x2 + 4)^2;
fx7 = (x1 + 4*x2 + 2)^2;
fx8 = (5*x1 + 4*x2 + 6)^2;

% 
% fz1 = chi(5,1)^2 ;
% fz2 = 2*chi(6,1)^2 - 3*chi(6,1) + 4;
% fz3 = 3*chi(7,1)^2 - 4*chi(7,1) + 3;
% fz4 = 4*chi(8,1)^2 - 5*chi(8,1) + 2;

%============== Gradient of fx============

% gd_fx11 = diff(fx1);
% gd_fx1x = gd_fx11(chi(1,1));
gd_fx1 = gradient(fx1);
gd_fx11 = subs(gd_fx1(1), [x1 x2], {chi(1,1),chi(2,1)});
gd_fx12 = subs(gd_fx1(2), [x1 x2], {chi(1,1),chi(2,1)});
gd_fx1 = double([gd_fx11 gd_fx12]');

gd_fx2 = gradient(fx2);
gd_fx21 =subs(gd_fx2(1), [x1 x2], {chi(3,1),chi(4,1)});
gd_fx22 =subs(gd_fx2(2), [x1 x2], {chi(3,1),chi(4,1)});
gd_fx2 = double([gd_fx21 gd_fx22]');

gd_fx3 = gradient(fx3);
gd_fx31 =subs(gd_fx3(1), [x1 x2], {chi(5,1),chi(6,1)});
gd_fx32 =subs(gd_fx3(2), [x1 x2], {chi(5,1),chi(6,1)});
gd_fx3 = double([gd_fx31 gd_fx32]');

gd_fx4 = gradient(fx4);
gd_fx41 =subs(gd_fx4(1), [x1 x2], {chi(7,1),chi(8,1)});
gd_fx42 =subs(gd_fx4(2), [x1 x2], {chi(7,1),chi(8,1)});
gd_fx4 = double([gd_fx41 gd_fx42]');

gd_fx5 = gradient(fx5);
gd_fx51 =subs(gd_fx5(1), [x1 x2], {chi(9,1),chi(10,1)});
gd_fx52 =subs(gd_fx5(2), [x1 x2], {chi(9,1),chi(10,1)});
gd_fx5 = double([gd_fx51 gd_fx52]');

gd_fx6 = gradient(fx6);
gd_fx61 =subs(gd_fx6(1), [x1 x2], {chi(11,1),chi(12,1)});
gd_fx62 =subs(gd_fx6(2), [x1 x2], {chi(11,1),chi(12,1)});
gd_fx6 = double([gd_fx61 gd_fx62]');

gd_fx7 = gradient(fx7);
gd_fx71 =subs(gd_fx7(1), [x1 x2], {chi(13,1),chi(14,1)});
gd_fx72 =subs(gd_fx7(2), [x1 x2], {chi(13,1),chi(14,1)});
gd_fx7 = double([gd_fx71 gd_fx72]');

gd_fx8 = gradient(fx8);
gd_fx81 =subs(gd_fx8(1), [x1 x2], {chi(15,1),chi(16,1)});
gd_fx82 =subs(gd_fx8(2), [x1 x2], {chi(15,1),chi(16,1)});
gd_fx8 = double([gd_fx81 gd_fx82]');




GDx = [gd_fx1; gd_fx2 ; gd_fx3 ; gd_fx4 ; gd_fx5 ; gd_fx6 ; gd_fx7; gd_fx8];

%============== Gradient of fz============
% gd_fz1 = 2*chi(5,1);
% gd_fz2 = 4*chi(6,1) - 3;
% gd_fz3 = 6*chi(7,1) - 4;
% gd_fz4 = 8*chi(8,1) - 5;
% 
% GDz = [gd_fz1 gd_fz2 gd_fz3 gd_fz4]';

%======== Laplacian =============
n = 8;
s1 = [1 1 2 2 2 3 3 4 4 4 4 5 5 5 5 6 6 7 7 8];
s2 = [2 5 1 4 3 2 4 3 2 5 8 1 4 6 7 5 7 5 6 4];
G = digraph(s1,s2);

in_deg = [];
for kk = 1 : n
    dd = find(s2 == kk);
    card_dd = length(dd);
    in_deg = [in_deg card_dd];
end
Adj = adjacency(G);
D = diag(in_deg);
L = D - Adj'; 

Fa = -eye(16,16);
Ba = [1 3 2 1 1 3 2 1;
      2 -1 3 2 1 3 2 1;
      -2 1 1 2 1 3 2 1;
      2 1 -4 3 1 3 2 1;
      2 1 0 0 0 0 0 0;
      2 0 0 0 0 0 0 0;
      2 1 0 0 1 0 0 0;
      0 1 0 0 1 0 0 0];
  Ba = kron(Ba, eye(2,2));

% %%%%%%%%%%%%%%%%%%%%%%%%%% Iqbal's virtual network%%%%%%%%%%%%%%%%%%
  alpha = 1/(1+t);
  
 L_bar = kron(L, eye(2,2));
 x_dot = - (L_bar *  chi(1:16,:) ) + (beta * L_bar * chi(17:32,:)) - (alpha *  GDx) +  (L_bar *  chi(33:48,:)); 
 z_dot =  -beta* ( L_bar * chi(1:16,:)) - ( L_bar * chi(17:32,:)) +  (L_bar *  chi(33:48,:))  ; %%%%%%% with (alpha *  GDz) LD is okay

 d_dot = Fa * chi(33:48,:) + (Ba * chi(1:16,:)) ;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  chi_dot = [x_dot ; z_dot ;  d_dot ];
end