

%% Initialization
clear;
close all
clc
% Environment parameters
c = 1500;       % speed of sound
f = 200;        % frequency
lambda = c/f;   % wavelength
set=[30 40 50 60 70];
iter = 5;
% SNR = 20;    % Signal-to-noise ratio

NMSE_set1 = zeros(length(set),iter);
NMSE_set2 = zeros(length(set),iter);
NMSE_set3 = zeros(length(set),iter);
NMSE_set4 = zeros(length(set),iter);
NMSE_set5 = zeros(length(set),iter);

jj = 1;
for m= set

    for MC = 1:iter

        rng(MC)
        n=100;                                          % signal dimension
       % m= 35;
        cp = 1;% number of measurements
        K=25;                                           % total number of nonzero coefficients
        L=3;                                            % number of nonzero blocks
        % rng(100)
        SNR =20;    % Signal-to-noise ratio
        % siga2=1e-12;
        Nsnapshot = 100

        % rng(5)
        % generate the block-sparse signal
        X=zeros(n,Nsnapshot);
        r=abs(randn(L,1)); r=r+1; r=round(r*K/sum(r));
        r(L)=K-sum(r(1:L-1));                           % number of non-zero coefficients in each block
        g=round(r*n/K);
        g(L)=n-sum(g(1:L-1));
        g_cum=cumsum(g);

        for i=1:L
            % generate i-th block
            seg = (1/sqrt(2)*(randn(r(i), Nsnapshot) +1j*randn(r(i),Nsnapshot)));              % generate the non-zero block
            % seg = 1*(1/sqrt(2)*(rand(r(i), Nsnapshot) +1i*rand(r(i),Nsnapshot)));              % generate the non-zero block
            %  seg = exp(1i*2*pi*randn(r(i), Nsnapshot))+1;              % generate the non-zero block
            % seg = randn(r(i), Nsnapshot);

            cp = 0.5;
            cp2 = 0;
            R = eye(r(i)) + diag(cp*ones(r(i)-1,1),1)+diag(cp*ones(r(i)-1,1),-1)+diag(cp2*ones(r(i)-2,1),2)+diag(cp2*ones(r(i)-2,1),-2);
            seg=sqrtm(R)*(seg);
            loc=randperm(g(i)-r(i));        % the starting position of non-zero block
            x_tmp=zeros(g(i), Nsnapshot);
            x_tmp(loc(1):loc(1)-1+r(i),:)= seg;
            X(g_cum(i)-g(i)+1:g_cum(i), :)=x_tmp;
        end

        % generate the measurement matrix
        Phi=randn(m,n);
        %A=Phi./(ones(m,1)*sqrt(sum(Phi.^2)));
        A=Phi;
        % noiseless measurements
        % X(X~=0) = 1;
        % cp=0.8;

        % X = X.*exp(1i*2*pi*randn(size(X)));
        % %






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

        % x_new=mu_new;
        % subplot(1,5,1)
        % nmse=norm(mean(abs(x_new),2)-mean(abs(X),2))^2/org_x
        %
        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % title(['EM SBL NMSE v1: ', num2str(nmse)])

        %% Revoery via PC-SBL
        eta=1;
        sigma2 =1;
        x_pcsbla=MPCSBL(Ysignal,A,sigma2,eta);
        nmse2(MC)=norm(mean(abs(x_pcsbla),2)-mean(abs(X),2))^2/org_x;

        % x_new=mu_new;
        % subplot(1,5,2)
        % nmse=norm(mean(abs(x_new),2)-mean(abs(X),2))^2/org_x
        %
        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % title(['EM PCSBL NMSE v1: ', num2str(nmse)])
        %% Revoery via PC-SBL
        eta=0;
        sigma2 =1;
        x_csbl0=MPCSBL_alternative(Ysignal,A,sigma2,eta);
        % x_new=mu_new;

        % org_x = norm(mean(abs(X),2))^2;
        nmse3(MC)=norm(mean(abs(x_csbl0),2)-mean(abs(X),2))^2/org_x;

        %% Revoery via PC-SBL
        eta=1;
        sigma2 =1;
        x_csbl=MPCSBL_alternative(Ysignal,A,sigma2,eta);
        % x_new=mu_new;

        % org_x = norm(mean(abs(X),2))^2;
        nmse4(MC)=norm(mean(abs(x_csbl),2)-mean(abs(X),2))^2/org_x;

        % subplot(1,6,1)
        %
        %
        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % legend('Ground Truth','Reconstructed');
        % title(['EM SBL NMSE: ', num2str(nmse)])


        %%

           %% Revoery via SBL
        eta=1;
        sigma2 =1;
        x_new=MSBL_correlated(y,A,sigma2,eta);
        % x_new=mu_new;

        nmse5(MC)=norm(mean(abs(x_fpcsbl),2)-mean(abs(X),2))^2/org_x;

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

    % nmse4(nmse4>nmse3) = 0;
    NMSE_set1(jj,:) = nmse1;
    NMSE_set2(jj,:) = nmse2;
    NMSE_set3(jj,:) = nmse3;
    NMSE_set4(jj,:) = nmse4;
    NMSE_set5(jj,:) = nmse5;

    jj = jj + 1;
end


figure

data = (1/100)*set;

plot(data,mean(NMSE_set1,2),'LineWidth',3)
hold on
plot(data,mean(NMSE_set2,2),'LineWidth',3)
hold on
plot(data,mean(NMSE_set3,2),'LineWidth',3)
hold on
plot(data,mean(NMSE_set4,2),'LineWidth',3)
hold on
plot(data,mean(NMSE_set5,2),'LineWidth',3)

ylabel('NMSE')
xlabel('m/n')
legend('EM SBL', 'EM CSBL','EM PCSBL','FP SBL', 'FP PCSBL' )
grid on
title(['K = ', num2str(K),  ', L = ', num2str(L), ', SNR = ', num2str(SNR)])
title(['Uncorrelated sources: K = ', num2str(K),  ', L = ', num2str(L), ', SNR = ', num2str(SNR)])


xlim([data(1) data(end)])