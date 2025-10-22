
%% Initialization
clear;
close all
clc
% Environment parameters
c = 1500;       % speed of sound
f = 200;        % frequency
lambda = c/f;   % wavelength
set=[20];
iter = 1;
% SNR = 20;    % Signal-to-noise ratio

NMSE_set1 = zeros(length(set),iter);
NMSE_set2 = zeros(length(set),iter);
NMSE_set3 = zeros(length(set),iter);
NMSE_set4 = zeros(length(set),iter);
NMSE_set5 = zeros(length(set),iter);

jj = 1;
for SNR= set

    for MC = 1

        % rng(1)
        s=randi([1,50])
        Nsnapshot = 1;
        Nsensor = 100;
        %SNR = 100;    % Signal-to-noise ratio
        sigma2=1e-12;
        % Nsensor =75;               % number of sensors
        d = 1/2*lambda;             % intersensor spacing
        q = (0:1:(Nsensor-1))';     % sensor numbering
        xq = (q-(Nsensor-1)/2)*d;   % sensor locations
        options = SBLSet();
        options.beta = 0;
        n=261;                                          % signal dimension
        m=Nsensor;                                           % number of measurements
        K=20;                                           % total number of nonzero coefficients
        L=3;
        % SNR =20;
        % number of nonzero blocks
        % generate the block-sparse signal
        x=zeros(n,1);
        r=abs(randn(L,1)); r=r+1; r=round(r*K/sum(r));
        r(L)=K-sum(r(1:L-1));                           % number of non-zero coefficients in each block
        g=round(r*n/K);
        g(L)=n-sum(g(1:L-1));
        %     rng(1)
        % sensor configuration structure
        Sensor_s.Nsensor = Nsensor;
        Sensor_s.lambda = lambda;
        Sensor_s.d = d;
        Sensor_s.q = q;
        Sensor_s.xq = xq;

        % SBL options structure for version 3

        options.convergence.error = 10^(-3);


        % total number of snapshots

        % number of random repetitions to calculate the average performance
        Nsim = 1;

        % range of angle space
        thetalim = [-90 90];

        theta_separation = 0.5;

        % Bearing grid
        theta = (thetalim(1):theta_separation:thetalim(2))';
        Ntheta = length(theta);

        % Design/steering matrix
        sin_theta = sind(theta);
        A = exp(-1i*2*pi/lambda*xq*sin_theta.')/sqrt(Nsensor);
        % A =randn(m,n+100);
        % to set Nsource parameter
        options.Nsource = K;

        %         rng(is, 'twister');
        g_cum=cumsum(g);

        for i=1:L
            % generate i-th block
             seg = 1/sqrt(2)*(rand(r(i), 1) +1j*rand(r(i),1));              % generate the non-zero block
            % seg = randn(r(i), 1);              % generate the non-zero block

             seg = 1/sqrt(2)*(randn(r(i), 1) +1j*randn(r(i),1));              % generate the non-zero block

            % seg = exp(1j*2*pi*randn(r(i), 1) );
            cp = 0.5;
            cp2 = 0;
            R = eye(r(i)) + diag(cp*ones(r(i)-1,1),1)+diag(cp*ones(r(i)-1,1),-1)+diag(cp2*ones(r(i)-2,1),2)+diag(cp2*ones(r(i)-2,1),-2);
            seg=sqrtm(R)*seg;

            loc=randperm(g(i)-r(i));        % the starting position of non-zero block
            x_tmp=zeros(g(i), 1);
            x_tmp(loc(1):loc(1)-1+r(i))= seg;
            x(g_cum(i)-g(i)+1:g_cum(i), 1)=x_tmp;
        end
        loc(s)
        cp =1;


        x = [zeros(50,1);x;zeros(50,1)];
        % x(x~=0) = 1;
        X = repmat(x,[1 Nsnapshot]);
        % X = X.*exp(1i*2*pi*randn(size(X)));

        %         X = X + cp* [zeros(1,Nsnapshot); X(1:end-1,:)]+cp* [X(2:end,:); zeros(1,Nsnapshot)];
        % cp = 0.8;
        % cp2 = 0;
        % R = eye(n) + diag(cp*ones(n-1,1),1)+diag(cp*ones(n-1,1),-1)+diag(cp2*ones(n-2,1),2)+diag(cp2*ones(n-2,1),-2);
        % X=sqrtm(R)*X;

        measure=A*X;

        % Observation noise, stdnoise = std(measure)*10^(-SNR/20);
        % stdnoise=sqrt(sigma2);
        % noise=randn(m,1)*stdnoise;

        % add noise to the signals
        rnl = 10^(-SNR/20)*norm(x);
        nwhite = complex(randn(Nsensor,Nsnapshot),randn(Nsensor,Nsnapshot))/sqrt(2*Nsensor);
        noise = nwhite * rnl;	% error vector

        % Noisy measurements
        Ysignal=measure+noise;


        org_x = norm(mean(abs(X),2))^2;

        %% Revoery via PC-SBL
        eta=0;
        sigma2 =1;
        x_ma=MPCSBL(Ysignal,A,sigma2,eta);
        % x_new=mu_new;

        org_x = norm(mean(abs(X),2))^2;

        nmse1(MC)=norm(mean(abs(x_ma),2)-mean(abs(X),2))^2/org_x;
        subplot(5,1,1)
        
        plot(mean(abs(X),2),'*')
        hold on
        stem(mean(abs(x_ma),2),'r')
        legend('Ground Truth','Reconstructed');
        title(['EM SBL NMSE: ', num2str(mean(nmse1))])


        %% Revoery via PC-SBL
        eta=1;
        sigma2 =1;
        tic
        x_pcsbla=MPCSBL(Ysignal,A,sigma2,eta);
        toc
        nmse2(MC)=norm(mean(abs(x_pcsbla),2)-mean(abs(X),2))^2/org_x;
        subplot(5,1,3)
        
        plot(mean(abs(X),2),'*')
        hold on
        stem(mean(abs(x_pcsbla),2),'r')
        % legend('Ground Truth','Reconstructed');
        title(['EM PCSBL NMSE: ', num2str(mean(nmse2))])
  
        %% Revoery via PC-SBL
        eta=1;
        sigma2 =1;
        x_csbl=MPCSBL_alternative(Ysignal,A,sigma2,eta);
        % x_new=mu_new;

        % org_x = norm(mean(abs(X),2))^2;
        nmse3(MC)=norm(mean(abs(x_csbl),2)-mean(abs(X),2))^2/org_x;
        subplot(5,1,2)

        plot(mean(abs(X),2),'*')
        hold on
        stem(mean(abs(x_csbl),2),'r')
        title(['EM CSBL NMSE: ', num2str(mean(nmse3))])
        % subplot(1,6,1)
        %
        %
        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % legend('Ground Truth','Reconstructed');
        % title(['EM SBL NMSE: ', num2str(nmse)])


        %% Revoery via PC-SBL
        options = SBLSet();
        options.Nsource = K;
        options.convergence.error   = 10^(-2);

        [x_fpsbl,report]=SBL( A , Ysignal, options );
        x_fpsbl = reshape(x_fpsbl,size(X))
        subplot(5,1,4)

        % x_new=mu_new;
        % subplot(1,6,2)
        nmse4(MC)=norm(mean(abs(x_fpsbl),2)-mean(abs(X),2))^2/org_x;

        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % title(['EM PCSBL NMSE: ', num2str(nmse)])
        plot(mean(abs(X),2),'*')
        hold on
        stem(mean(abs(x_fpsbl),2),'r')
        % figure
        % plot(report.results.error)
        title(['FP SBL NMSE: ', num2str(mean(nmse4))])

%%
        options.beta=0.9;
        tic
        [x_fpcsbl,report]=SBL_PC( A , Ysignal, options );
        toc
        x_fpcsbl = reshape(x_fpcsbl,size(X));
        % x_new=mu_new;
        subplot(5,1,5)
        nmse5(MC)=norm(mean(abs(x_fpcsbl),2)-mean(abs(X),2))^2/org_x;
        
        plot(mean(abs(X),2),'*')
        hold on
        stem(mean(abs(x_fpcsbl),2),'r')
        title(['FP PCSBL NMSE: ', num2str(mean(nmse5))])
        % x_new(abs(x_new)<0.4) = 0;
        % subplot(1,6,5)
        %
        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % hold on
        % plot(abs(diag(report.results.final_iteration.gamma)),'g')
        % hold on
        % plot(abs(diag(report.results.final_iteration.gamma,1)),'p')
        % legend('ground truth','reconstructed x','diagonal \gamma','subdiagonal \gamma')
        % title(['FP TriSBL NMSE: ', num2str(nmse)])

    end

   
end


