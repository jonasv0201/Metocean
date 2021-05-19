clear;close;clc;
load('s8');
%% Analyse distribution of most probable largest crest height in each storm.
%Keep only struct 
%clearvars  -except s8
L = s8.stormMaxCrest.likelihood;
crest = s8.stormMaxCrest.crestHeight;
crest = sortrows(crest);
clearvars -except crest s8
paramMom=methodOfMoments(crest);
paramMLE=MLE(crest);
paramWBL2 = weibull2(crest);

s8.dots = {};
s8.dots.mom=paramMom;
s8.dots.mle=paramMLE;
s8.dots.wlb2=paramWBL2;
%% Functions
function params = weibull2(data)
gamma = min(data)-0.1;
params = wblfit(data-gamma);
params(3) = gamma;
end
function params = methodOfMoments(data)

[mean, var, skew] = calculateMoments(data);
save('tempMoments','mean','var','skew');
fun = @wbl3Params;
%params = [alfa, beta, lambda, mean, var, skew]
params = fsolve(fun,[1,1,5]);
params(4:6)=[mean,var,skew];
end
function F = wbl3Params(x)
load('tempMoments','mean','var','skew');

G1 =@(b) gamma( 1 + (1./b) );
G2 =@(b) gamma( 1 + (2./b) );
G3 =@(b) gamma( 1 + (3./b) );
meanFunc = @(a,b,c) c + a.*G1(b);
varFunc = @(a,b) (a.^2) .* ( G2(b) - (G1(b).^2) );
skewFunc = @(b) ( G3(b) - 3.*G1(b).*G2(b) + 2.*(G1(b).^3) ) ./ ( G2(b) - G1(b).^2 ).^(1.5);
%Set of functions F = 0
F(1) = mean - meanFunc(x(1),x(2),x(3));
F(2) = var - varFunc(x(1),x(2));
F(3) = skew - skewFunc(x(2));
end
function [Mean, Var, Skew ]= calculateMoments(vec)
n = length(vec);
Mean = mean(vec);
Var = var(vec,0); 
skewfunc = @(x) (sum((x-Mean).^3)./n) ./ Var.^1.5;
Skew = skewfunc(vec);
end
function params = MLE(data,initalGuess)
if ~exist('initialGuess','var')
      initialGuess = [2.2 1.3 7.9];
end
wbl3 = @(x,a,b,c) (x>c).*(b/a).*(((x-c)/a).^(b-1)).*exp(-((x-c)/a).^b);

opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'FunValCheck','off');

params = mle(data,'pdf',wbl3,'start',initialGuess,...
     'Options',opt,'LowerBound',[0.1 0.1 5],'UpperBound',[10 10 min(data)]);
end