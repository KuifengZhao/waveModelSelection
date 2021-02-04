# waveModelSelection
Basic water wave theory
This function is to be used together with the waveTypeChart.fig and the dataset of criteria - waveZoneData.mat. In the chart, range of h/L is [10^(-2), 1] and a/L is [2 x 10^(-6), 0.1]. h is water depth, e.g. 1m; T is wave period, e.g., 2s; H is wave height = 2*a = 0.17m; % wave amplitude a = H/2;

waveRange is the zone in the redrawn chart and waveType is suggested wave model.

Eg 1, 

[waveRange, waveType]  = waveModelSelection('h', 1, 'T', 2, 'H',0.17);

waveType =

  1×1 cell array
    {'Stokes III'}
    
Eg 2,

waveH = [0.17 0.17 0.02 0.02];
waveh = [1 2 1 2];
waveT = [2 2 2 2];

[waveRange, waveType] = waveModelSelection('H', waveH, 'h', waveh, 'T', waveT)

waveType =

  4×1 cell array

    {'Stokes III'}
    {'Stokes II' }
    {'Stokes II' }
    {'Laminar'   }
