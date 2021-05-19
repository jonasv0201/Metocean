clear;close;clc;
load('s8');
%% Analyse distribution of most probable largest crest height in each storm.
%Keep only struct 
%clearvars  -except s8
crest = s8.stormMaxCrest.crestHeight;
MOM = s8.dots.mom(1:3);
WBL2 = s8.dots.wbl2;
MLE = s8.dots.mle;

x=7.5:0.01:17;
wbl3 = @(x,a,b,c) (x>c).*(b/a).*(((x-c)/a).^(b-1)).*exp(-((x-c)/a).^b);

grid on; hold on;
histogram(crest,'BinWidth',1,'lineStyle','none','FaceColor','green','Normalization','probability','FaceAlpha',0.2)
histogram(crest,'BinWidth',0.2,'lineStyle','none','FaceColor','red','Normalization','probability','FaceAlpha',0.8)

plot(x,wbl3(x,MOM(1),MOM(2),MOM(3)),'--','LineWidth',1.3);
%Almost identical to MLE
%plot(x,wbl3(x,WBL2(1),WBL2(2),WBL2(3)),'--','LineWidth',1.3); 
plot(x,wbl3(x,MLE(1),MLE(2),MLE(3)),'LineWidth',1.3,'Color','blue');
title('Most probable maximum crest height - Weibull 3 fit');
xlabel('Crest height [m]')
ylabel('PDF'); 
legend('Histogram, bin size 1.0','Histogram, bin size 0.2','MoM, Weibull3','MLE, Weibull3')