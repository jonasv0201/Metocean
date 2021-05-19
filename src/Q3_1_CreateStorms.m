clear;close;
load('windSeaData.mat','rawData')
% storms
% 1    2        3        4        5          6     7    8    9      10
% Hs  Tp_adj  Tp_old   Windsped  WindDir   Year  Month  Day  Hour  index
lim = 8;
j = 1;
%% Filter out storms above 8m
for i = 1:length(rawData)
if( rawData(i,1)>lim )
   states(j,1:9)  = rawData(i,:);
   states(j,10) = i; %add Index
    j = j+1; 
end
end
states( all(~states,2), : ) = []; %Delete all zero rows.


Z=zeros(281,12);
s8 = struct('hs',Z,'tp',Z,'t1',Z,'af',Z,'bf',Z,'s1',Z,'ur',Z,'k1',Z,'idx',Z,'t2',Z);

%% Create struct of all the storms. 
%S8.Hs(4,7) is Hs ofstorm number 4, step 7
num = 1;
step = 1;
s8.hs(num,step) = states(1,1);
s8.tp(num,step) = states(1,2);
s8.idx(num,step) = states(1,10);
for i = 2:length(states)
if((s8.idx(num,step)+1) == states(i,10))
    step = step + 1;
    s8.hs(num,step) = states(i,1);
    s8.tp(num,step) = states(i,2);
    s8.idx(num,step) = states(i,10);
else
    step = 1;
    num = num + 1;
    s8.hs(num,step) = states(i,1);
    s8.tp(num,step) = states(i,2);
    s8.idx(num,step) = states(i,10);
end
end
%Keep only struct 
clearvars  -except s8
%% Cycle trough each Hs&Tp pair and create Forristall numbers. DNVGL Rp205 3.5.10
step=0;
for i = 1:281
 while (s8.hs(i,step+1)) %While step is not zero
    step=step+1;
    Hs = s8.hs(i,step);
    Tp = s8.tp(i,step);
    [af,bf,k1,T1,s1,Ur,T2] = forristallParams(Hs,Tp);
    s8.af(i,step) = af;
    s8.bf(i,step) = bf;
    s8.k1(i,step) = k1;
    s8.t1(i,step) = T1;
    s8.s1(i,step) = s1;
    s8.ur(i,step) = Ur;
    s8.t2(i,step) = T2;
 end
    step=0;
end
%Keep only struct 
clearvars  -except s8
function [af,bf,k1,T1,s1,Ur,T2] = forristallParams(Hs,Tp)
% (8.70) in Metocean book / 3.5.10 In DNVGL RP205
D = 350;  % meters. Haltenbanken.
G = 9.81; 
funAf = @(s1,ur) 0.3536 + 0.2892.*s1 + 0.1060.*ur;
funBf = @(s1,ur) 2 - 2.1597.*s1 + 0.0968.*ur.^2;
funS1 = @(Hs, T1) 6.283.*Hs./(G.*T1.^2);
funUr = @(Hs, k1) Hs./((k1^2).*(D^3));
funT1 = @(Tp,peakShape) Tp.*(5+peakShape)./(6.8+peakShape);
funK1 = @(T1) 4*pi^2/(G.*T1^2);
funT2 = @(Tp,peakShape) Tp.*sqrt((5+peakShape)./(11+peakShape));

   T1 = funT1(Tp, pshape(Hs,Tp) );
   k1 = funK1(T1);
   s1 = funS1(Hs,T1);
   Ur = funUr(Hs,k1);
   bf = funBf(s1,Ur);
   af = funAf(s1,Ur);
   T2 = funT2(Tp, pshape(Hs,Tp) );
end
function gamma = pshape(Hs,Tp)
%DNVGL RP205 3.5.5.5
%Peak shape parameter, gamma
param = Tp/sqrt(Hs);
if (param <= 3.6)
    gamma = 5;
elseif (param<5)
   gamma = exp(5.75-1.15.*param);
else
    gamma = 1;
end %if
end %func