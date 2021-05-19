clear;close;clc;
load('s8');
%% 8.70 in metocean book
% Forristall Weibull  distribution of crest height
% CDF
n3h = @(T2) 10800/T2;
F_Crest = @(x,HsXaf,bf) (1 - exp(- ( x./HsXaf ).^bf ));
F_Max_C = @(x,HsXaf,bf,n3h) ( F_Crest(x,HsXaf,bf)).^n3h;

ddx_F_Crest = @(x,HsXaf,bf) (1./HsXaf^bf).*(exp(- ( x./HsXaf).^bf)).*bf.*x.^(bf-1) ;
ddx_F_Max =@(x,HsXaf,bf,n3h) n3h*(F_Crest(x,HsXaf,bf)).^(n3h-1).*ddx_F_Crest(x,HsXaf,bf);

step = 0;
x = 0.1:0.01:25;
first = find(x==6);
last = find(x==22);
figure();
hold on;
grid on;
PDF=zeros(281,length(x));
CDF_stormMaxCrest = ones(281,length(x));
for i = 1:281
 sumTerm = zeros(1,length(x));
 while (s8.hs(i,step+1)) %While step is not zero
    step=step+1;
    HsXaf = (s8.hs(i,step)*s8.af(i,step));
    bf=s8.bf(i,step);
    N = n3h(s8.t2(i,step));
    F_Max_Crest_Temp = F_Max_C(x, HsXaf, bf, N);
    CDF_stormMaxCrest(i,:) = CDF_stormMaxCrest(i,:).*F_Max_Crest_Temp;
    sumTerm = sumTerm + (ddx_F_Max(x,HsXaf,bf,N)./F_Max_Crest_Temp);
 end %While
 PDF(i,:) = CDF_stormMaxCrest(i,:) .* sumTerm;
 plot(x(first:last),CDF_stormMaxCrest(i,first:last));
 step=0;
end

%% Find Most probable Largest Crest Height for each storm
likelihood = zeros(1,281);
crestHeight = zeros(1,281);
CDFy= zeros(1,281);
idx = 0;
for i = 1:281
  idx = find(PDF(i,:)==max(PDF(i,:)));
  likelihood(i) = PDF(i,idx);
  crestHeight(i) = x(idx);
  CDFy(i) = CDF_stormMaxCrest(i,idx);
end
scatter(crestHeight,CDFy,'.','k');
%% Plot PDF
title('Storm maximum crest height')
xlabel('Crest height [m]')
ylabel('Cumulative distribution')
figure()
hold on;
grid on;
title('Storm maximum crest height')
xlabel('Crest height [m]')
ylabel('Probability density')
for i = 1:281
    plot(x(first:last),PDF(i,first:last));
end
scatter(crestHeight,likelihood,'.','k');

%% save Variables to S8
s8.stormMaxCrest = {};
s8.stormMaxCrest.PDF = PDF(:,first:last);
s8.stormMaxCrest.CDF = CDF_stormMaxCrest(:,first:last);
s8.stormMaxCrest.xdata = x(1,first:last);
s8.stormMaxCrest.likelihood = likelihood;
s8.stormMaxCrest.crestHeight = crestHeight;

%Keep only struct 
clearvars  -except s8

