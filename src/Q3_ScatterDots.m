clear;close;clc;
load('s8');
%% Scatter of most probable maximum Crest for each storm.
%Keep only struct 
%clearvars  -except s8
figure();
scatter(1:281,s8.stormMaxCrest.crestHeight,'.','k');
title('Most probable maximum crest height in each Storm');
ylabel('Crest height [m]');
xlabel('Storm number');
grid on;
%% Scatter Wave height vs max Hs.
figure();
peakHs = zeros(1,281);
for i=1:281
   peakHs(i) = max(s8.hs(i,:)); 
end
scatter(peakHs,2.*s8.stormMaxCrest.crestHeight,'.','k');
title('Most probable maximum wave height in each Storm');
ylabel('Wave height [m]');
grid on;
xlabel('Largest Hs in each storm [m]');