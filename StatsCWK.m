
%% ME3 STATS COURSEWORK
%AROHAN SUBRAMONIA
rng(01054062);
 
clear all
%% Reading the data from the csv file
AS10415 = readtable('as10415.csv');
DATA=table2cell(AS10415);
num=cell2mat(DATA(:,1:4));
num(:,5) = string(DATA(:,5)) == 'petrol'; %Fuel type logical (0 if Diesel, 1 if Petrol)
num(:,6) =  string(DATA(:,6)) == 'white';  %Colour logical (1 if white, 0 if other)

% Splitting the data up into four sections
l100=num(:,1);
mass=num(:,2);
t100=num(:,3);
disp=num(:,4);

%% Question 1
%% Exploratory Data Analysis
A_mean(1:4) = mean(num(:,1:4)); %Arithmetic mean
G_mean(1:4) = geomean(num(:,1:4)); %Geometric mean
median(1:4) = median(num(:,1:4)); %Median
T_mean(1:4) = mean(num(51:450,1:4)); %10% trimmed mean
A_std(1:4)  = std(num(:,1:4)); %Arithmetic standard deviation
top(1:4)    = sum(log(num(:,1:4)/G_mean).^2); %Geometric standard deviation
G_std(1:4)  = exp((top/499).^(0.5)); 

%% Histogram of Fuel Efficiency
histogram(num(:,1));
grid ON
fs=20; %FontSize
%title('Histogram of Fuel Efficiency per 100km','FontSize',18);
set(gca,'FontSize',fs)
xlabel('Fuel Efficiency (litres/100km)','FontSize',fs);
ylabel('Frequency','FontSize',fs);

%% Scatterplots
subplot(3,2,1), scatter(l100,mass,'.','k'), xlabel({'Fuel Efficiency';' (litres/100km)'},'FontSize',fs), ylabel('Vehicle Mass (kg)','FontSize',fs), grid MINOR
subplot(3,2,2), scatter(l100,t100,'.','k'), xlabel({'Fuel Efficiency';' (litres/100km)'},'FontSize',fs), ylabel({'Time taken to reach';' 100km/h (s)'},'FontSize',fs), grid MINOR
subplot(3,2,3), scatter(l100,disp,'.','k'), xlabel({'Fuel Efficiency';' (litres/100km)'},'FontSize',fs), ylabel('Engine Size (litres)','FontSize',fs), grid MINOR
subplot(3,2,4), scatter(mass,t100,'.','k'), xlabel('Vehicle Mass (kg)','FontSize',fs), ylabel({'Time taken to reach';' 100km/h (s)'},'FontSize',fs), grid MINOR
subplot(3,2,5), scatter(mass,disp,'.','k'), xlabel('Vehicle Mass (kg)','FontSize',fs), ylabel('Engine Size (litres)','FontSize',fs), grid MINOR
subplot(3,2,6), scatter(t100,disp,'.','k'), xlabel({'Time taken to reach';' 100km/h (s)'},'FontSize',fs), ylabel('Engine Size (litres)','FontSize',fs), grid MINOR

%% Boxplot drawing routine
fueltype=[{['Petrol'] ['Diesel']}; {['White'] ['Other']}]; %Box labels

XLabel=[{'Fuel Type'} {'Fuel Color'}];
YLabel=[{['Fuel Efficiency'];[' (litres/100km)']};{'Vehicle Mass';'(kg)'};{'Time taken to';'reach 100km/h (s)'};{'Engine size';'(litres)'}];
yl=[1,1,3,3,5,5,7,7];
var=[1,1,2,2,3,3,4,4];
fs=18;
b=[1,2,1,2,1,2,1,2];
for a=1:8
    subplot(4,2,a), boxplot(num(:,var(a)),num(:,b(a)+4),'labels',fueltype(b(a),:));
    xlabel(XLabel(b(a)),'FontSize',fs);
    ylabel(YLabel(yl(a):yl(a)+1),'FontSize',fs);   
    set(gca,'FontSize',fs)
end

%% Question 2
%% Standardizing and Organising data
Nl100=zscore(num(:,1));
Nmass=zscore(num(:,2));
Nt100=zscore(num(:,3));
Ndisp=zscore(num(:,4));
type = num(:,5);
colour=num(:,6);
%Creating a table of standardized variables
Nvarnames(1:6)=[{'Nl100'} {'Nmass'} {'Nt100'} {'Ndisp'} {'type'} {'colour'}];
Nnum=table(Nl100,Nmass,Nt100,Ndisp,num(:,5),num(:,6),'VariableNames',Nvarnames);

%% Linear Models - single variables without interaction

%Linear model of l100 vs mass
LMmass1=fitlm(mass,l100,'linear');
plot(LMmass1);
set(gca,'FontSize',fs)
ylabel('Fuel Efficiency (litres/100km)','FontSize',fs);
xlabel('Vehicle Mass (kg)','FontSize',fs);
title(' ')

%Quadratic model of l100 vs mass
LMmass2=fitlm(mass,l100,'purequadratic');
subplot(), plot(LMmass2);
set(gca,'FontSize',fs)
ylabel('Fuel Efficiency (litres/100km)','FontSize',fs);
xlabel('Vehicle Mass (kg)','FontSize',fs);

%Linear model of l100 vs t100
LMtime1=fitlm(t100,l100,'linear');
plot(LMtime1);
set(gca,'FontSize',fs)
ylabel('Fuel Efficiency (litres/100km)','FontSize',fs);
xlabel('Time taken to reach 100km/h (s)','FontSize',fs);

%Quadratic model of l100 vs t100
LMtime2=fitlm(t100,l100,'purequadratic');
plot(LMtime2);
set(gca,'FontSize',fs)
ylabel('Fuel Efficiency (litres/100km)','FontSize',fs);
xlabel('Time taken to reach 100km/h (s)','FontSize',fs);

%Linear model of l100 vs disp
LMdisp1=fitlm(disp,l100,'linear');
plot(LMdisp1);
set(gca,'FontSize',fs)
ylabel('Fuel Efficiency (litres/100km)','FontSize',fs);
xlabel('Engine Size (litres)','FontSize',fs);

%Quadratic model of l100 vs disp
LMdisp2=fitlm(disp,l100,'purequadratic');
plot(LMdisp2);
set(gca,'FontSize',fs)
ylabel('Fuel Efficiency (litres/100km)','FontSize',fs);
xlabel('Engine Size (litres)','FontSize',fs);

%% Parameters for single variable models without interaction (not standardized)
%Parameters for l100 vs mass - linear
Mass_linear=fitlm(Nmass,Nl100,'linear');
RSm1=Mass_linear.Rsquared.Ordinary;
MSEm1=Mass_linear.MSE;
AICm1=Mass_linear.ModelCriterion.AIC;

%Parameters for l100 vs mass - quadratic
Mass_quad=fitlm(Nmass,Nl100,'purequadratic');
RSm2=Mass_quad.Rsquared.Ordinary;
MSEm2=Mass_quad.MSE;
AICm2=Mass_quad.ModelCriterion.AIC;

%Parameters for l100 vs time - linear
Time_linear=fitlm(Nt100,Nl100,'linear');
RSt1=Time_linear.Rsquared.Ordinary;
MSEt1=Time_linear.MSE;
AICt1=Time_linear.ModelCriterion.AIC;

%Parameters for l100 vs time - quadratic
Time_quad=fitlm(Nt100,Nl100,'purequadratic');
RSt2=Time_quad.Rsquared.Ordinary;
MSEt2=Time_quad.MSE;
AICt2=Time_quad.ModelCriterion.AIC;

%Parameters for l100 vs disp - linear
Disp_linear=fitlm(Ndisp,Nl100,'linear');
RSd1=Disp_linear.Rsquared.Ordinary;
MSEd1=Disp_linear.MSE;
AICd1=Disp_linear.ModelCriterion.AIC;

%Parameters for l100 vs disp - quadratic
Disp_quad=fitlm(Ndisp,Nl100,'purequadratic');
RSd2=Disp_quad.Rsquared.Ordinary;
MSEd2=Disp_quad.MSE;
AICd2=Disp_quad.ModelCriterion.AIC;

RES_RS=[RSm1;RSm2;RSt1;RSt2;RSd1;RSd2];
RES_MSE=[MSEm1;MSEm2;MSEt1;MSEt2;MSEd1;MSEd2];
RES_AIC=[AICm1;AICm2;AICt1;AICt2;AICd1;AICd2];

%Final table of single variable parameters
modelname={'mass (linear)','mass (quad)','time (linear)','time (quad)','engine size (linear)','engine size (quad)'}';
estnames={'Type' 'Rsquared' 'MSE' 'AIC'};
statistics=table(modelname,RES_RS,RES_MSE,RES_AIC,'VariableNames',estnames);

%% Linear Models - multivariable, some with interactions

%Linear multivariate model, no interactions
Linear_Model = fitlm(Nnum,'linear','PredictorVars',{'Nmass','Nt100','Ndisp','type','colour'},'ResponseVar','Nl100');
lintab=table({'Linear Multivariate'},Linear_Model.Rsquared.Ordinary,Linear_Model.MSE,Linear_Model.ModelCriterion.AIC,'VariableNames',estnames);

%Quadratic multivariate model, no interactions
Quad_Model = fitlm(Nnum,'purequadratic','PredictorVars',{'Nmass','Nt100','Ndisp','type','colour'},'ResponseVar','Nl100');
quadtab=table({'Pure Quadratic Multivariate'},Quad_Model.Rsquared.Ordinary,Quad_Model.MSE,Quad_Model.ModelCriterion.AIC,'VariableNames',estnames);

%Linear multivariate model, with interactions
Inter_Model1 = fitlm(Nnum,'interactions','PredictorVars',{'Nmass','Nt100','Ndisp','type','colour'},'ResponseVar','Nl100');

%Quadratic multivariate model, with interactions
Inter_Model2 = fitlm(Nnum,'quadratic','PredictorVars',{'Nmass','Nt100','Ndisp','type','colour'},'ResponseVar','Nl100');

inter1tab=table({'Linear Interactions'},Inter_Model1.Rsquared.Ordinary,Inter_Model1.MSE,Inter_Model1.ModelCriterion.AIC,'VariableNames',estnames);
inter2tab=table({'Quadratic Interactions'},Inter_Model2.Rsquared.Ordinary,Inter_Model2.MSE,Inter_Model2.ModelCriterion.AIC,'VariableNames',estnames);

%Final table of all models
statistics(7:10,1:4)=[lintab;quadtab;inter1tab;inter2tab];

%% Refining the best performing (quadratic) model

LM2=Inter_Model2;
%Removing the least significant coefficients
for a=1:(0.5*length(LM2.Coefficients.Estimate))
    [X,I]=min(abs(LM2.Coefficients.Estimate));
    if I ~= 1
        LM2=removeTerms(LM2,char(LM2.CoefficientNames(I)));
    end
end

for a=1:10
    LM2=step(LM2);
end

%% Adding terms in to finalise model
mdl_eqn2 = ('Nl100~1+Nmass+Nt100+type+Nmass:Nt100+Nmass:type+Nt100:type');
LM2=fitlm(Nnum,mdl_eqn2);

%% Question 3
%% Residual and Fitted Plot
res=LM2.Residuals.Raw;
fit=LM2.Fitted;

scatter(fit,res,10,'k','MarkerFaceColor','k','MarkerEdgeColor','k');
grid MINOR;
set(gca,'FontSize',fs);
xlabel('Fitted Values','FontSize',fs);
ylabel('Residuals','FontSize',fs);

%% QQ Plot of residuals
qqplot(res);
grid MINOR;
title('Q-Q Plot of Residuals from model','FontSize',fs);
set(gca,'FontSize',fs);

%% Coefficient Plot
coeff=LM2.Coefficients.Estimate;

plot(coeff,'k.-','MarkerSize',15);
grid MINOR;
set(gca,'FontSize',fs);
xlabel('Coefficients index','FontSize',fs);
ylabel('Coefficients values','FontSize',fs);

%% Question 4
%% Bootstrapping
mint=min(t100);
maxt=max(t100);
time=[mint:0.1:maxt];
res1=LM2.Residuals.Raw;
l100hat=LM2.Fitted;
res=l100-l100hat;

%Using IMstar, predict 
for j=1:10
    res_star = randsample(res,length(res),'true');
    l100star=LM2.Fitted+res_star;
    
    % Same model as Q2 fitted to Nl100star
    mdl_eqnstar = ('l100star~1+mass+t100+type+mass:t100+mass:type+t100:type');
    numstar=table(l100star,mass,t100,disp,num(:,5),num(:,6),'VariableNames',{'l100star' 'mass' 't100' 'disp' 'type' 'colour'});
    LM2star=fitlm(numstar,mdl_eqnstar);

    for i=1:length(time)
        l100starhat(i)=LM2star.Coefficients.Estimate(1)+...
            LM2star.Coefficients.Estimate(2).*mean(mass)+...
            LM2star.Coefficients.Estimate(3).*time(i)+...
            LM2star.Coefficients.Estimate(4).*min(type)+...
            LM2star.Coefficients.Estimate(5).*(mean(mass).*time(i))+...
            LM2star.Coefficients.Estimate(6).*(mean(mass).*min(type))+...
            LM2star.Coefficients.Estimate(5).*(time(i).*min(type));
    end

    l100stardata(:,j)=l100starhat';
end

l100starmean=mean(l100stardata');
lower = quantile(l100stardata,0.025,2);
upper = quantile(l100stardata,0.975,2);

scatter(time,l100starmean,5,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold;
plot(time,lower,'k.-',time,upper,'k.-');
grid MINOR
set(gca,'FontSize',fs)
xlabel('Time Array (s)','FontSize',fs);
ylabel({'Predicted Fuel Efficiencies';'(litres/100km)'},'FontSize',fs);

%% Using a separate model to predict confidence intervals

xnew=zeros(length(time),5);
xnew(:,1)=mean(mass);
xnew(:,2)=time;
xnew(:,3)=mean(disp);

[ypred,yci]=predict(LM2star,xnew);

plot(time,ypred,'r.-',time,yci,'r.-');
grid MINOR
set(gca,'FontSize',fs)
xlabel('Time Array (s)','FontSize',fs);
ylabel({'YPred Fuel Efficiencies';'(litres/100km)'},'FontSize',fs);

%% Question 5
%% Reading the data from the csv file
testdata = readtable('testdata.csv');
TDATA=table2cell(AS10415);
tnum=cell2mat(TDATA(:,1:4));
tnum(:,5) = string(TDATA(:,5)) == 'petrol'; %Fuel type logical (0 if Diesel, 1 if Petrol)
tnum(:,6) =  string(TDATA(:,6)) == 'white';  %Colour logical (1 if white, 0 if other)

% Splitting the data up into six sections
Tl100=tnum(:,1);
Tmass=tnum(:,2);
Tt100=tnum(:,3);
Tdisp=tnum(:,4);
Ttype=tnum(:,5);
Tcolour=tnum(:,6);

%% Using Question 2 model on Test data
Tvarnames(1:6)=[{'Tl100'} {'Tmass'} {'Tt100'} {'Tdisp'} {'Ttype'} {'Tcolour'}];
Tnum=table(Tl100,Tmass,Tt100,Tdisp,Ttype,Tcolour,'VariableNames',Tvarnames);

mdl_eqnT = ('Tl100~1+Tmass+Tt100+Ttype+Tmass:Tt100+Tmass:Ttype+Tt100:Ttype');
LMT=fitlm(Tnum,mdl_eqnT);

%Quadratic Interactions model on test data
Test_Inter_Quad_Model=fitlm(Tnum,'quadratic','PredictorVars',{'Tmass' 'Ttype' 'Tt100' 'Tdisp' 'Ttype' 'Tcolour'},'ResponseVar','Tl100');
