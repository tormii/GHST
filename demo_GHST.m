git  
%仿真信号TFR不同方法对比
clc; clear all; close all;    
    N = 2000;
    fs = 200;
    t = (0:N-1)/fs;
    f = (0:N/2)*fs/N;
  
%% 信号生成
%模式1
    a1 = 0;
    b1 = 1;
    c1 = 0;
    d1 = 0;
    e1 = 0;
    A_f1 = exp(0.005*f);
    GD_t1 = b1+c1*f+d1*f.^2+e1*f.^3;
    Phi_f1 = a1+b1*f+c1*f.^2/2+d1*f.^3/3+e1*f.^4/4;    
    X1 = A_f1.*exp(-1i*2*pi*Phi_f1);
    X1(end) = -A_f1(end);
    Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
    y1 = ifft(Y1);
%模式2
    a2 = 0;
    b2 = 2;
    c2 = 1/50;
    d2 = 0;
    e2 = 0;
    A_f2 = exp(0.008*f);
    GD_t2 = b2+c2*f+d2*f.^2+e2*f.^3;
    Phi_f2 = a2+b2*f+c2*f.^2/2+d2*f.^3/3+e2*f.^4/4;    
    X2 = A_f2.*exp(-1i*2*pi*Phi_f2);
    X2(end) = -A_f2(end);%
    Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
    y2 = ifft(Y2);
%模式3
    a3 = 0;
    b3 = 31/10;
    c3 = 7/50;
    d3 = -1/1000;
    e3 =0;
    A_f3 = exp(0.01*f);
    GD_t3 = b3+c3*f+d3*f.^2+e3*f.^3;
    Phi_f3 = a3+b3*f+c3*f.^2/2+d3*f.^3/3+e3*f.^4/4;    
    X3 = A_f3.*exp(-1i*2*pi*Phi_f3);
    X3(end) = -A_f3(end);%
    Y3 = [X3  conj(fliplr(X3(2:end-1)))];    
    y3 = ifft(Y3);
%加噪声
figure;
plot(GD_t1,f,GD_t2,f,GD_t3,f);
    axis tight; 
    xlabel('Time (t)','FontSize',10);
    ylabel('Frequency (Hz)','FontSize',10);
    xlim([0 10]);
y = y3;
figure;plot(t,real(y));
    axis tight; 
    xlabel('Time (t)','FontSize',10);
    ylabel('Amplitude','FontSize',10);
%     y = y2;
%      y = awgn(y,noise,'measured');         %为信号加载白噪声

%% 输入参数
    noise = 0;             %信号高斯噪声分贝，值越小噪声越大
    fs = 200;              %信号采样率
    N = 2000;               %信号点数
    %GHST计算阶数选择
    Mode = '3Ord';        %(1Ord，2Ord)
    %窗选择参数
    WindowOpt = struct('type','gauss','s',0.12);%0.05
    %频率选择参数
    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%        y = awgn(y,noise,'measured');         %为信号加载白噪声
    %%  STFT 
    [Wx,t1,f1,~] = stft(y, fs, WindowOpt, Parameter, 'normal');
    figure;
    imagesc(t1,f1,abs(Wx));axis xy;
    title('STFT','FontSize',14);
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Frequency (Hz)','FontSize',10);
%      hold on;plot([7.8 8.3],[58 58]);
%     hold on;plot([7.8 8.3],[78 78]);
%     hold on;plot([7.8 7.8],[58 78]);
%     hold on;plot([8.3 8.3],[58 78]);
%     r = renyi(abs(Wx),t,f',3)
%% GHST
    [hst,t,f,xMean,~] = GHST(y , fs,  WindowOpt, Parameter, Mode);
    figure;
    imagesc(t,f,abs(hst));axis xy;
    title('GHST','FontSize',14);
    axis tight; xlabel('Time (t)','FontSize',10);
    ylabel('Frequency (Hz)','FontSize',10);
%     xlim([0.9 1.1]);ylim([80 85]);
%     xlim([3.6 3.7]);ylim([80 85]);
%     xlim([7.8 7.9]);ylim([80 85]);
%     hold on;plot([7.8 8.3],[58 58]);
%     hold on;plot([7.8 8.3],[78 78]);
%     hold on;plot([7.8 7.8],[58 78]);
%     hold on;plot([8.3 8.3],[58 78]);
%     r = renyi(abs(hst),t,f',3)

%%  重构 
   singleside=15;
   direction=3;
    Real_GD=GD_t3;
    Ex = ExtractTSST_new( hst, Real_GD, fs,singleside,direction);
%     figure;
%    imagesc(t,f,abs(Ex));axis xy;
%    title('GHSTEX','FontSize',14);
%    axis tight; xlabel('Time (t)','FontSize',10);
%    ylabel('Frequency(Hz)','FontSize',10);
   
    b = itsst(Ex,fs, xMean);
    figure
    plot(t,y3,t,b);
  
    
    title('Reconstruction','FontSize',14);
    axis tight; xlabel('Time (t)','FontSize',10);
    ylabel('AMP','FontSize',10);
    legend('Original','New');
   
    