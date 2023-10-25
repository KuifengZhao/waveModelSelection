clc,clear
% this code works in Matlab version R2020a, earlier version not guarantted.
kh = (10.^(-3:0.02:0.8))'; % arbitrary h/L range in typical Le Mehaute Chart
% hOverL = kh/2/pi;
kH = 10.^(-5:0.02:-0.23)*2; % arbitrary H/L range in typical Le Mehaute Chart
theta = 0:0.01:2*pi;
[HOverL, hOverL] = meshgrid(kH/2/pi, kh/2/pi);
axRange = [0.01 1 1e-5 0.3];
URLim = 26; % Stokes wave applicable in the range with Ursell <=26
HOverL(HOverL>(URLim.*hOverL.^3))= nan;
%% given kH and kh, solve for ka
syms kasym
B31 = (3+8*(tanh(kh)).^2-9*(tanh(kh)).^4    )/16./(tanh(kh)).^4;
B33 = (27-9*(tanh(kh)).^2+9*(tanh(kh)).^4-3*(tanh(kh)).^6)/64./(tanh(kh)).^6;
B51 = (121*cosh(2*kh).^5 + 263*cosh(2*kh).^4 + 376*cosh(2*kh).^3 - 1999*cosh(2*kh).^2 ...
    +2509*cosh(2*kh) - 1108)./(192*(cosh(2*kh) - 1).^5);
B53 = (57*cosh(2*kh).^7 + 204*cosh(2*kh).^6 - 53*cosh(2*kh).^5 ...
    -782*cosh(2*kh).^4 - 741*cosh(2*kh).^3 - 52*cosh(2*kh).^2 ...
    + 371*cosh(2*kh) + 186)*9./((3*cosh(2*kh) + 2).*(cosh(2*kh) - 1).^6*128); 
B55 = (300*cosh(2*kh).^8 + 1579*cosh(2*kh).^7 + 3176*cosh(2*kh).^6 + 2949*cosh(2*kh).^5 ...
    + 1188*cosh(2*kh).^4 + 675*cosh(2*kh).^3 + 1326*cosh(2*kh).^2 + 827*cosh(2*kh) ...
    + 130)*5./((cosh(2*kh) - 1).^6.*(12*cosh(2*kh).^2 + 11*cosh(2*kh) + 2)*384);
%%
kHfun = 2*pi*HOverL -kasym.*(2+2*(B31+B33).*kasym.^2+2*(B51+B53+B55).*kasym.^4) ==0;

for i = 1:length(kh)
    for j = 1:length(kH) 
        if ismissing(HOverL(i,j))
            ka(i,j) =nan;
        else
            kasol = double(vpasolve(kHfun(i,j), [HOverL(i,j)/3*2*pi inf]));
            ka(i,j) = kasol(1);
        end
    end
end

%% plot each order magnitude components of free surface 
sigma = tanh(kh);
alpha1 = cosh(2.*kh);

Hovera = 2*(1+ (B31+B31).*ka.^2+(B51+B55+B53).*ka.^4);

for i = 1:length(theta)
% first and second order
eta1(:,:,i) = ones(length(kh),length(ka))*cos(theta(i)); 
eta2(:,:,i) = ka./4.*(3-sigma.^2)./sigma.^3 .* cos(2.*theta(i));
% our 3rd order

eta3(:,:,i) = ka.^2.*(B31.*cos(theta(i))+B33.*cos(3.*theta(i)));
% our derived 4th order equation
sigma1 =24.*(3.*cosh(2.*kh)+2).*(cosh(2.*kh)-1).^4.*sinh(2.*kh);
eta4(:,:,i) = ka.^3./sigma1.*((60.*alpha1.^6 ...
+232.*alpha1.^5-118.*alpha1.^4-989.*alpha1.^3-607.*alpha1.^2 ...
+352.*alpha1+260).*cos(2.*theta(i))+(24.*alpha1.^6 ...
+116 .*alpha1.^5+214.*alpha1.^4+188.*alpha1.^3+133.*alpha1.^2 ...
+101.*alpha1 +34).*cos(4.*theta(i)));
% our derived 5th order equation
eta5(:,:,i) = ka.^4.*B51.*cos(theta(i)) ...
+ka.^4.*B53.*cos(3.*theta(i)) ...
+ka.^4.*B55.*cos(5.*theta(i));
end

eta = eta1+eta2+eta3+eta4+eta5;
etaRatio2 = eta2(:,:,1)./eta1(:,:,1); %Definition of R2 in the updated manuscript submitted, under review
etaRatio3 = eta3(:,:,1)./(eta1(:,:,1)+eta2(:,:,1));
etaRatio4 = eta4(:,:,1)./(eta1(:,:,1)+eta2(:,:,1)+eta3(:,:,1));
etaRatio5 = eta5(:,:,1)./(eta1(:,:,1)+eta2(:,:,1)+eta3(:,:,1)+eta4(:,:,1));

% etaRatio2 = max(abs(eta2),[],3)./max(abs(eta1),[],3); %superseded by the results evaluated at crest Line64. 
% etaRatio3 = max(abs(eta3),[],3)./(max(abs(eta1),[],3)+max(abs(eta2),[],3));
% etaRatio4 = max(abs(eta4),[],3)./(max(abs(eta1),[],3)+max(abs(eta2),[],3)+max(abs(eta3),[],3));
% etaRatio5 = max(abs(eta5),[],3)./(max(abs(eta1),[],3)+max(abs(eta2),[],3)+max(abs(eta3),[],3)+max(abs(eta4),[],3));

%% A2 over A1
% HoverL = Hovera.*aOverL;
figure
[C2,h2] = contourf(hOverL, HOverL, etaRatio2, [0.01 0.05 0.1 0.2],'showText', 'on');
colormap(jet(6))
caxis([0 0.3])
xlabel('h/L'), ylabel('H/L')
set(gca, 'XScale', 'log', 'YScale', 'log')
axRange = [0.01 1 2e-6 0.2];
axis(axRange)

C2_cid = find(C2(1,:)==0.01); %the 1% contour line
A2overA1_1perc =[C2(1,2:1+C2(2,1));C2(2,2:1+C2(2,1))]'

%% A3 over A3prec
figure
[C3,h3] = contourf(hOverL, HOverL,etaRatio3, [0.01 0.05 0.1 0.2],'showText', 'on');
colormap(jet(6))
caxis([0 0.3])

set(gca, 'XScale', 'log', 'YScale', 'log')
% axRange = [0.01 1 2e-6 0.2];
axis(axRange)
C3_cid = find(C3(1,:)==0.01); %the 1% contour line
A3overA3p_1perc =  [C3(1,2:1+C3(2,1));C3(2,2:1+C3(2,1))]';

%% A4 over A4prec
figure
[C4,h4] = contourf(hOverL, HOverL, etaRatio4, [0.01 0.05 0.1 0.2],'showText', 'on');
colormap(jet(6))
caxis([0 0.3])

set(gca, 'XScale', 'log', 'YScale', 'log')
axRange = [0.01 1 2e-6 0.2];
axis(axRange)
C4_cid = find(C4(1,:)==0.01); %the 1% contour line
% C4 - column 207:end is for ratio = 0.05 and higher, 155-207 is truncated as a very
% narrow wedge in contour plot
A4overA4p_1perc = [C4(1,84:155); ...
    C4(2,84:155)]';
%% A5 over A5prec
figure
[C5,h5] = contourf(hOverL, HOverL, etaRatio5, [0.01 0.05 0.1 0.2],'showText', 'on');
colormap(jet(6))
caxis([0 0.3])

set(gca, 'XScale', 'log', 'YScale', 'log')
axRange = [0.01 1 2e-6 0.2];
axis(axRange)
C5_cid = find(C5(1,:)==0.01); %the 1% contour line
A5overA5p_1perc =[C5(1,2:1+C5(2,1));C5(2,2:1+C5(2,1))]';

%% combined plots
figure
hold on,
plot(A2overA1_1perc(:,1), A2overA1_1perc(:,2),'k')
plot(A3overA3p_1perc(:,1), A3overA3p_1perc(:,2),'k') 
plot(A4overA4p_1perc(:,1), A4overA4p_1perc(:,2),'k') 
plot(A5overA5p_1perc(:,1), A5overA5p_1perc(:,2),'k') 
set(gca, 'XScale', 'log', 'YScale', 'log')
% colorbar
ax = gca;
ax.TickLength = [0.03, 0.03];
axis on,
set(gcf, 'Color', 'w')
set(gca, 'Layer', 'top')
% Breaking and ursell number
hold on,
% Breaking Lines by Fenton
lambdaOverh = 2*pi./kh;
HoverLBreak = kh.*(0.141063.*lambdaOverh+0.0095721.*lambdaOverh.^2+0.0077829.*lambdaOverh.^3)...
    ./(1+0.0788340.*lambdaOverh+0.0317567.*lambdaOverh.^2+0.0093407.*lambdaOverh.^3)/2/pi;
plot(hOverL, HoverLBreak, 'b-.', 'linewidth', 2)
axRange = [0.01 1 1e-5 0.3];
axis(gca, axRange);
set(gca, 'XScale', 'log', 'YScale', 'log')
 
% creteria for Ursell number
Ur10 = 10*(kh/2/pi).^3; % line of Ursell No = 1;
Ur10(Ur10>HoverLBreak) = nan;
plot(kh/2/pi, Ur10, 'b', 'linewidth', 2)
Ur1 = 1*(kh/2/pi).^3; % line of Ursell No = 1;
Ur1(Ur1>HoverLBreak) = nan;
plot(kh/2/pi, Ur1, 'b--', 'linewidth', 2)
Ur26 = 26*(kh/2/pi).^3; % line of Ursell No = 1;
Ur26(Ur26>HoverLBreak) = nan;
plot(kh/2/pi, Ur26, 'b--', 'linewidth', 2)

% colorbar
ax = gca;
axis(gca, axRange);
ax.TickLength = [0.03, 0.03];
%% plot m curves based on Fenton's cnoidal wave solutions
% Case for given H, h, m
h = 10; % mean water depth, the line for other depth was tested to be identical with 
H = 0.01:0.001:7.8;

% m =1-4e-16;
m = 0.96;
[Km, Em] = ellipke(m);
em = Em./Km;
Hoverh = H/h;
% Fentons (1999) - Eqn B.7 in his work on cnoidal wave
L = (4*Km./sqrt(3*Hoverh).*(1+Hoverh.*(5/8-3/2.*em)+...
    Hoverh.^2.*(-21/128+1/16.*em+3/8.*em.^2)+Hoverh.^3.*(20127/179200-...
    409/6400.*em+7/64.*em.^2+1/16.*em.^3)+...
    Hoverh.^4.*(-1575087/28672000+1086367/1792000.*em-2679/25600.*em.^2+...
    13/128.*em.^3+3/128.*em.^4))).*h;
plot(h./L, H./L, 'r:','linewidth',2)

m =1-4e-8;
[Km, Em] = ellipke(m);
em = Em./Km;
Hoverh = H/h;
% Fentons (1999) - Eqn B.7 in his work on cnoidal wave
L = (4*Km./sqrt(3*Hoverh).*(1+Hoverh.*(5/8-3/2.*em)+...
    Hoverh.^2.*(-21/128+1/16.*em+3/8.*em.^2)+Hoverh.^3.*(20127/179200-...
    409/6400.*em+7/64.*em.^2+1/16.*em.^3)+...
    Hoverh.^4.*(-1575087/28672000+1086367/1792000.*em-2679/25600.*em.^2+...
    13/128.*em.^3+3/128.*em.^4))).*h;
plot(h./L, H./L, 'r','linewidth',2)
%% plot numerical solutions region as analytical solution becomes problematic
solLim = kh/2/pi*0.4; solLim (kh/2/pi >0.20004) = nan;
solLim (kh/2/pi <0.024) = nan;
hold on, plot(kh/2/pi, solLim, 'k-.','linewidth',2)
vLine(0.5);
vLine(0.05);
tightfig
%% annotation

xt = [ 1.2e-2, 2.1e-2, 2.7e-2];
yt = [6e-5, 6e-5,1.4e-5];
annoRot = {'$U_r = 26$','$U_r = 10, m \approx 0.5$', '$U_r = 1$'};
text(xt, yt, annoRot, 'Rotation', 45, 'Color', 'b', 'Interpreter', 'Latex');
xlabel('$h/L$', 'Interpreter', 'Latex'), ylabel('H/L', 'Interpreter', 'Latex')

