function [waveParamZone, waveModel] = waveModelSelection(varargin)

% This function is to be used together with the waveModelChart.fig and the
% dataset of criteria - waveModelData.mat. 
% In the chart, range of h/L is [10^(-2), 1] and a/L is [2 x 10^(-6), 0.1]
% h is water depth, e.g. 1m;
% T is wave period, e.g., 2s; 
% H is wave height = 2*a = 0.17m; % wave amplitude a = H/2;
% waveType is suggested wave theory type
% Eg 1, 
% [waveParamZone, waveModel] = waveModelSelection('h', 1, 'T', 2, 'H',0.17);
% waveParamZone =
% 
%   1×1 cell array
%     {'Stokes III'}
% Eg 2,
% waveH = [0.17 0.17 0.02 0.02];
% waveh = [1 2 1 2];
% waveT = [2 2 2 2];
% [waveParamZone, waveModel] = waveModelSelection('H', waveH, 'h', waveh, 'T', waveT)
% waveParamZone =
% 
%   4×1 cell array
% 
%     {'Stokes III'}
%     {'Stokes II' }
%     {'Stokes II' }
%     {'Laminar'   }

h = 1; %water depth, m
H = 0.17; %wave height, m
T = 2; %wave period, s
names = varargin(1:2:end);
values = varargin(2:2:end);
for kIn = 1:numel(names)
    switch names{kIn}
        case 'H'
            H = values{kIn};
        case 'h'
            h = values{kIn};
        case 'T'
            T = values{kIn};
    end
end

% coefficients for boundary
C1 = 10/2; C2 = 0.01; C3 = 1/8; C4 = 1/2;
% Guo (2002)'s formula about dispersion relationship is used in this code
g = 9.81;
hgT2 = h./g./(T.^2);
hoverL = 2*pi*hgT2.*(1-exp(-(4*pi^2.*hgT2).^(5/4))).^(-2/5);
L = h./hoverL; 
aoverL = H/2./L;
open('waveModelChart.fig')
hold on, plot(hoverL, aoverL, 'ks','markerface','w')

% check range of data, which polygon is the points
load('waveModelData.mat')
% hoverLMax = 1;
% aoverLMax = 0.07;
waveModel = cell(length(h),1); waveParamZone = cell(length(h),1);
% make polygon for linear range
waveLinear = [A2overA1_1perc; ...
    [A2overA1_1perc(1,1),A2overA1_1perc(end,2)];...
    A2overA1_1perc(1,:)];
[in1,on1] = inpolygon(hoverL,aoverL, waveLinear(:,1),waveLinear(:,2));
waveParamZone(or(in1,on1)) = {'A'}; %1;
waveModel(or(in1,on1)) = {'Linear fully dispersive equation (LFDE) or Stokes 1st order equation'}; %1;
% make polygon for Stokes II range
waveStokes2 = [A2overA1_1perc; ...
    flip(A3overA3p_1perc);...
    A2overA1_1perc(1,:)];
[in2,on2] = inpolygon(hoverL,aoverL, waveStokes2(:,1),waveStokes2(:,2));
waveParamZone(or(in2,on2)) = {'B'}; %1;
waveModel(or(in2,on2)) = {'2nd order'}; %2;
% make polygon for Stokes III range
waveStokes3 = [A3overA3p_1perc; ...
    flip(A5overA5p_1perc);...
    A3overA3p_1perc(1,:)];
[in3,on3] = inpolygon(hoverL,aoverL, waveStokes3(:,1),waveStokes3(:,2));
waveParamZone(or(in3,on3)) = {'C'}; %1;
waveModel(or(in3,on3)) = {'3rd order'}; %3;
% make polygon for Stokes V range
kh =(0.001:0.01:6.3)';
hOverL = kh/2/pi;
lambdaOverh = 2*pi./kh;
kaFenton = kh/2.*(0.141063.*lambdaOverh+0.0095721.*lambdaOverh.^2+0.0077829.*lambdaOverh.^3)...
    ./(1+0.0788340.*lambdaOverh+0.0317567.*lambdaOverh.^2+0.0093407.*lambdaOverh.^3);

limMiche = [A5overA5p_1perc(:,1), 0.07*tanh(A5overA5p_1perc(:,1)*2*pi)];
waveStokes5 = [A5overA5p_1perc; ...
    flip(limMiche);...
    A5overA5p_1perc(1,:)];
[in5,on5] = inpolygon(hoverL,aoverL, waveStokes5(:,1),waveStokes5(:,2));
waveParamZone(or(in5,on5)) = {'D'}; %1;
waveModel(or(in5,on5)) = {'5th order'}; %5;
% check range of Long wave;
% Ur10 = C1*hOverL.^3; % line of Ursell No = 10;
% Ur10limLow = InterX([hOverL'; aOverh], [hOverL'; Ur10']);

aOverh = C2*hOverL;
UrUpper = C1*hOverL.^3; % line of Ursell No = 10;
UrLower = C4*hOverL.^3; % line of Ursell No = 10;
RangeE1 = find(hOverL >1/20);
RangeF1 = find(hOverL >C3);

UrUpperlimLow = InterX([hOverL'; aOverh'], [hOverL'; UrUpper']);
UrUpplim = InterX([hOverL'; kaFenton'/2/pi], [hOverL'; UrUpper']);
RangeI1 = find(hOverL > UrUpperlimLow(1,1));
RangeI2 = find(hOverL > UrUpplim(1,1));
% E linear non dispersive shallow water wave
LNDSWW = [[hOverL(1:RangeE1(1)-1), UrLower(1:RangeE1(1)-1)];...
    [hOverL([RangeE1(1)-1, 1]), UrLower([1,1])]];
[inLNDSWW,onLNDSWW] = inpolygon(hoverL,aoverL, LNDSWW(:,1),LNDSWW(:,2));
waveParamZone(or(inLNDSWW,onLNDSWW)) = {'E'}; %1;
waveModel(or(inLNDSWW,onLNDSWW)) = {'Linear non-dispersive equation (LNDE) or linear shallow water wave (LSWE)'}; %7;

% F linear weakly dispersive shallow water wave
LWDSWW = [[hOverL(1:RangeF1(1)-1), UrLower(1:RangeF1(1)-1)];...
    [hOverL([RangeF1(1)-1, 1]), UrLower([1,1])]];
[inLWDSWW,onLWDSWW] = inpolygon(hoverL,aoverL, LWDSWW(:,1),LWDSWW(:,2));
waveParamZone(or(inLWDSWW,onLWDSWW)) = {'F'}; %1;
waveModel(or(inLWDSWW,onLWDSWW)) = {'Linear weakly dispersive shallow water wave (LWDE)'}; %7;

Boussi = [[hOverL(1:RangeF1(1)-1), UrLower(1:RangeF1(1)-1)];...
    [hOverL(RangeF1(1)-1), aOverh(RangeF1(1)-1)];...
    flip([hOverL(RangeI1(1):RangeF1(1)-1), aOverh(RangeI1(1):RangeF1(1)-1)]);...
    flip([hOverL(1:RangeI1(1)-1), UrUpper(1:RangeI1(1)-1)])];
[inBoussi,onBoussi] = inpolygon(hoverL,aoverL, Boussi(:,1),Boussi(:,2));
waveParamZone(or(inBoussi,onBoussi)) = {'G'}; %1;
waveModel(or(inBoussi,onBoussi)) = {'Classic Boussinesq approximation'}; %7;


LnonDSWW = [[hOverL(1:RangeI1(1)-1), UrUpper(1:RangeI1(1)-1)];...
    flip([hOverL(1:RangeI1(1)-1), aOverh(1:RangeI1(1)-1)]);...
    [hOverL(1), UrUpper(1)]];
[inLnonDSWW,onLnonDSWW] = inpolygon(hoverL,aoverL, LnonDSWW(:,1),LnonDSWW(:,2));
waveParamZone(or(inLnonDSWW,onLnonDSWW)) = {'H'}; %1;
waveModel(or(inLnonDSWW,onLnonDSWW)) = {'Linear non-dispersive equation (LNDE) or linear shallow water wave (LSWE)'}; %7;

NLSWW = [[hOverL(1:RangeI1(1)-1), aOverh(1:RangeI1(1)-1)]; ...
    [hOverL(RangeI1(1):RangeI2(1)-1), UrUpper(RangeI1(1):RangeI2(1)-1)];...
    flip([hOverL(1:RangeI2(1)), kaFenton(1:RangeI2(1))/2/pi]);[hOverL(1), UrUpper(1)]];
[inNLSWW,onNLSWW] = inpolygon(hoverL,aoverL, NLSWW(:,1),NLSWW(:,2));
waveParamZone(or(inNLSWW,onNLSWW)) = {'I'}; %1;
waveModel(or(inNLSWW,onNLSWW)) = {'Non-linear shallow water equation (NLSWE)'}; %7;
% check breaking region
yaxisRange = ylim; 

waveBreak = [[hOverL, kaFenton/2/pi];...
    [hOverL(end), yaxisRange(2)];[hOverL(1), yaxisRange(2)];...
    [hOverL(1), kaFenton(1)/2/pi]];
[inBreak,onBreak] = inpolygon(hoverL,aoverL, waveBreak(:,1),waveBreak(:,2));
waveParamZone(or(inBreak,onBreak)) = {'J'}; %8;
waveModel(or(inBreak,onBreak)) = {'Breaking'}; %8;

