This function is to be used together with the waveModelChart.fig and the dataset of criteria - waveModelData.mat. In the chart, range of h/L is [10^(-2), 1] and a/L is [2 x 10^(-6), 0.1]. h is water depth, e.g. 1m; T is wave period, e.g., 2s; H is wave height = 2*a = 0.17m; % wave amplitude a = H/2;

waveRange is the zone in the redrawn chart and waveType is suggested wave model. A detailed explanation of the wave model selection process can be found from the working manuscript 'A guide for selection of water wave models'.

Eg 1, 

[waveParamZone, waveModel]  = waveModelSelection('h', 1, 'T', 2, 'H',0.17);

waveModel =

  1×1 cell array
    {'Stokes III'}
    
Eg 2,

waveH = [0.17 0.17 0.02 0.02];
waveh = [1 2 1 2];
waveT = [2 2 2 2];

[waveParamZone, waveModel] = waveModelSelection('H', waveH, 'h', waveh, 'T', waveT)

waveModel =

  4×1 cell array

    {'Stokes III'}
    {'Stokes II' }
    {'Stokes II' }
    {'Laminar'   }
