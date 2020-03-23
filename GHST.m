function [Tx,t,f,xMean,GD] = GHST(x , fs,  WindowOpt, Parameter, Mode)
%	------------------- Generalized Horizontal Synchrosqueezing Transform -------------------- 
% Authors: Xiaotong Tu and Fucai Li
% email:tormiier@gmail.com,tormii@sjtu.edu.cn;
% https://www.researchgate.net/profile/Xiaotong_Tu2
%Input:
%       
%       x:imput signal
%       fs:sampling frequency/采样频率（Hz）
%       WindowOpt:window function/窗函数选择参数
%           WindowOpt.s：(0.01) 窗函数初始尺度
%           WindowOpt.f0：(0) 窗函数初始中心频率
%           WindowOpt.type：(gauss) 窗函数类型
%       Parameter:频率选择参数
%           Parameter.L：(200) 频率划分个数
%           Parameter.fmin：(最小分辨率) 分析最小频率
%           Parameter.fmax：(奈奎斯特频率) 分析最大频率
%       Mode:(1Ord，2Ord,3Ord)
%Output:
%       Tx:TFR
%       t:time
%       f:frequency/压缩后的频率（Hz）
%       GroupDelay:group delay;
%---------------------------------------------------------------------------------
% When using this code, please cite our paper:
% Xiaotong Tu,Qi Zhang, Zhoujie He, Yue Hu, Saqlain Abbas and Fucai Li, Generalized Horizontal Synchrosqueezing Transform: Algorithm and Applications, IEEE Transactions on Industrial Electronics
% Author: Xiaotong Tu（Mar.,2020）
%---------------------------------------------------------------------------------
%% 预处理信号
    N = length(x);
%% 参数赋值
    s = WindowOpt.s; type = WindowOpt.type;
    L = Parameter.L; fmin = Parameter.fmin; fmax = Parameter.fmax;
    gamma =sqrt(eps); 
    %% SST计算
    %STFT计算
    [Wx00,t,f,xMean] = stft(x, fs, WindowOpt, Parameter, 'normal');
     %群延迟
    if strcmp(Mode, '1Ord')
        WindowOpt.type = '1ord(w)_gauss';
        [Wx01,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        GD = real(Wx01./Wx00/(1i));
        for ptr = 1:N
            GD(:,ptr) = GD(:,ptr) + t(ptr);
        end
        GD( abs(Wx00) < gamma ) = Inf;
    elseif strcmp(Mode, '2Ord')
        WindowOpt.type = '1ord(w)_gauss'; 
       [Wx01,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
       WindowOpt.type = 'w1*1ord(w)_gauss'; 
       [Wx11,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
       WindowOpt.type = 'w1*gauss'; 
       [Wx10,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
       WindowOpt.type = 'w2*gauss'; 
       [Wx20,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
 
         A1=-Wx11.*Wx00+Wx01.*Wx10-Wx00.*Wx00;%第二次代码
         A2=Wx20.*Wx00-Wx10.*Wx10;
         x21=Wx10./Wx00;
         y1= -Wx01./Wx00/(1i);
         for ptr = 1:N
            y1(:,ptr) = y1(:,ptr) - t(ptr);
         end
           y2=A1./(1i*A2);
           temp2=y2;
%            save('temp2.mat','temp2');
%            temp2=y2.*x21;
           GD=real(-y1+y2.*x21);
           GD( abs(Wx00) < gamma ) = Inf;
           GD( abs(A2) < gamma ) = Inf;
    elseif strcmp(Mode, '3Ord')
       WindowOpt.type = '1ord(w)_gauss'; 
       [Wx01,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
       WindowOpt.type = 'w1*1ord(w)_gauss'; 
       [Wx11,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = 'w2*1ord(w)_gauss'; 
       [Wx21,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = 'w1*gauss'; 
       [Wx10,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
       WindowOpt.type = 'w2*gauss'; 
       [Wx20,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
       WindowOpt.type = 'w3*gauss'; 
       [Wx30,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
       WindowOpt.type = 'w4*gauss'; 
       [Wx40,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');

%第二次代码

     
          A1=-Wx11.*Wx00+Wx01.*Wx10-Wx00.*Wx00;%第二次代码 
          dA1=1i*(-Wx21.*Wx00+Wx01.*Wx20-2*Wx00.*Wx10);
          A2=Wx20.*Wx00-Wx10.*Wx10;
          dA2=1i*(Wx30.*Wx00-Wx10.*Wx20);
          B1=Wx30.*Wx00-Wx20.*Wx10;
          dB1=1i*(Wx40.*Wx00-Wx20.*Wx20);
          x21=Wx10./Wx00;
          x31=Wx20./Wx00;
          x32=B1./A2;
          y1= -Wx01./Wx00/(1i);
         for ptr = 1:N
            y1(:,ptr) = y1(:,ptr) - t(ptr);
         end
           y2=A1./(1i*A2);
           temp3=y2;
%             save('temp3.mat','temp3');
           y3=(dA1.*A2-A1.*dA2)./(1i*(dB1.*A2-B1.*dA2));
%             GD=real(-y1+y2.*x21);% N=2
%             temp3=x21.*(y2-x32.*y3)+x31.*y3;
            GD=real(-y1+x21.*(y2-x32.*y3)+x31.*y3);
           GD( abs(Wx00) < gamma ) = inf;
           GD( abs(A2) < gamma ) = inf;
            GD( abs(dB1.*A2-B1.*dA2) < gamma ) = inf;
    else
        error('Unknown SST Mode: %s', Mode);
   
    end
    %频率差分计算
    dt = 1/fs;
    %窗时域生成
    [gf,~] = windowf(s,type);
    %计算g(0)
    g0 = gf(0);
    g0 = conj(g0);
    if(g0 == 0)
        error('window must be non-zero and continuous at 0 !');
    end
    %同步压缩
    Wx00(isinf(GD)) = 0;
    Tx = zeros(L,N);
    
    for prt=1:L
        for b=1:N
           
            m = min(max(1 + round((GD(prt,b)-0)/dt),1),N);
            Tx(prt, m) = Tx(prt, m) + Wx00(prt, b)*dt;
        end
    end
    
    Tx = Tx / g0;
end