  clc;clear all;close all;
  
%% Coherence
Fs=200;
y1=T4_segment3;
y2=F4_segment3;
[S1,F,T]=spectrogram(y1,Fs*0.5,Fs*0.45,4000,Fs);
% 412 59 49 512

%-----Normalize
[aa,bb]=size(S1);
Smax=max(S1);
for k=1:bb
    for kk=1:aa
        S1(kk,k)=S1(kk,k)/Smax(k);
    end
end
delta=mean(abs(S1(1:8,:)));%delta 0.1-3.9Hz18
theta=mean(abs(S1(9:16,:)));%theta 4-7.9Hz
alpha=mean(abs(S1(17:24,:)));%alpha 8-11.9Hz
beta=mean(abs(S1(25:60,:)));%beta 12-29.9Hz
slowgamma=mean(abs(S1(51:101,:)));%slowgamma 25-50Hz61101
gamma=mean(abs(S1(61:241,:)));%gamma 30-120Hz61101
All=mean(abs(S1(1:242,:)));
delta_1=smoothdata(delta,'gaussian',20);
theta_1=smoothdata(theta,'gaussian',20);
alpha_1=smoothdata(alpha,'gaussian',20);
beta_1=smoothdata(beta,'gaussian',20);
slowgamma_1=smoothdata(slowgamma,'gaussian',20);
gamma_1=smoothdata(gamma,'gaussian',20);
All_1=smoothdata(All,'gaussian',20);
%
%Fs=20000;
[S2,F,T]=spectrogram(y2,Fs*0.5,Fs*0.45,4000,Fs);
%-----Normalize
[aa2,bb2]=size(S2);
Smax2=max(S2);
for k2=1:bb2
    for kk2=1:aa2
        S2(kk2,k2)=S2(kk2,k2)/Smax2(k2);
    end
end
delta=mean(abs(S2(1:8,:)));%delta 0.1-3.9Hz
theta=mean(abs(S2(9:16,:)));%theta 4-7.9Hz
alpha=mean(abs(S2(17:24,:)));%alpha 8-11.9Hz
beta=mean(abs(S2(25:60,:)));%beta 12-29.9Hz
slowgamma=mean(abs(S2(51:101,:)));%slowgamma 25-50Hz
gamma=mean(abs(S2(61:241,:)));%gamma 30-120Hz
All=mean(abs(S2(1:242,:)));
delta_2=smoothdata(delta,'gaussian',20);
theta_2=smoothdata(theta,'gaussian',20);
alpha_2=smoothdata(alpha,'gaussian',20);
beta_2=smoothdata(beta,'gaussian',20);
slowgamma_2=smoothdata(slowgamma,'gaussian',20);
gamma_2=smoothdata(gamma,'gaussian',20);
All_2=smoothdata(All,'gaussian',20);
%%
[C121,lag121]= xcorr(delta_1,delta_2,'coeff');
 rG1_delta=max(C121);
   [M21,I21] = max(C121);
   t21 = lag121(I21);
subplot(2,1,1)
plot(lag121,C121,[t21 t21],[-0.5 1],'r:')
 text(t21+100,0.5,['Lag: ' int2str(t21) 'point'])
 text(t21+100,0,['r= '  num2str(rG1_delta)])
 ylabel('Clus-delta')
 axis tight
 title('Cross-Correlations')
 box off
%%
[C121,lag121]= xcorr(slowgamma_1,slowgamma_2,'coeff');
 rG1_slowgamma=max(C121);
   [M21,I21] = max(C121);
   t21 = lag121(I21);
subplot(2,1,1)
plot(lag121,C121,[t21 t21],[-0.5 1],'r:')
 text(t21+100,0.5,['Lag: ' int2str(t21) 'point'])
 text(t21+100,0,['r= '  num2str(rG1_slowgamma)])
 ylabel('Clus-slowgamma')
 axis tight
 title('Cross-Correlations')
 box off
%%
[C121,lag121]= xcorr(gamma_1,gamma_2,'coeff');
 rG1_gamma=max(C121);
   [M21,I21] = max(C121);
   t21 = lag121(I21);
subplot(2,1,1)
plot(lag121,C121,[t21 t21],[-0.5 1],'r:')
 text(t21+100,0.5,['Lag: ' int2str(t21) 'point'])
 text(t21+100,0,['r= '  num2str(rG1_gamma)])
 ylabel('Clus-gamma')
 axis tight
 title('Cross-Correlations')
 box off
 %%
[C121,lag121]= xcorr(All_1,All_2,'coeff');
 rG1_All=max(C121);
   [M21,I21] = max                                              (C121);
   t21 = lag121(I21);
subplot(2,1,1)
plot(lag121,C121,[t21 t21],[-0.5 1],'r:')
 text(t21+100,0.5,['Lag: ' int2str(t21) 'point'])
 text(t21+100,0,['r= '  num2str(rG1_All)])
 ylabel('Clus-All')
 title('Cross-Correlations')
 box off
%%
% %  yyaxis left
plot(T,gamma_1,'k','LineWidth',1)
hold on 
%  yyaxis right
plot(T,gamma_2,'r','LineWidth',1)
box off
%%
% %  yyaxis left
plot(T,delta_1,'k','LineWidth',1)
hold on 
%  yyaxis right
plot(T,delta_2,'r','LineWidth',1)
box off