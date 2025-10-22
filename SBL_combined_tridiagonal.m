function [ x_post, report] = SBL_combined_tridiagonal_mv( A , Y, options )
%
% function [ gamma , report ] = SBL_v3p1( A , Y, options )
% The idea behind SBL is to find a diagonal replica 'covariance' Gamma.
% Minimizing (YY^T / AGA^T + penality) should lead to the correct
% replica selection (up to a bogus scale factor/amplitude).
%
% Attention: If Y is single snapshot (and single frequency), it needs to be
% a row vector (the code makes a 2nd snapshot with repmat).
%
% Inputs
%
% A - Multiple frequency augmented dictionary <f , n , m>
%     f: number of frequencies
%     n: number of sensors
%     m: number of replicas
%   Note: if f==1, A = < n , m >
%
%
% Y - Multiple snapshot multiple frequency observations <f , n , L>
%     f: number of frequencies
%     n: number of sensors
%     L: number of snapshots
%   Note: if f==1, Y = < n , L >
%
% options - see SBLset.m
%
%
% Outputs
%
% gamma <m , 1> - vector containing source power
%                 1: surfaces found by minimum error norm
%
% report - various report options
%

%% check function

if ismatrix(A) % SBL needs frequency dimension
    B(1,:,:) = A;
    A        = B;
end

% number of frequencies
Nfreq     = size(A,1);

% single frequency single snapshot
if ismatrix(Y)
    % either 1 freq or 1 snapshot
    if  Nfreq == 1
        
        if size(Y,2) == 1 % single snapshot
            Y=Y.'; %squeeze2
        else
            Y = permute(Y,[ 3 1 2 ]); % works
        end
        
    end
    
end

%%

options.SBL_v = '3.12';

%% slicing

Nsource = options.Nsource;

if options.tic == 1
    tic
end

%% Initialize variables

% number of sensors
Nsensor   = size(A,2);
% number of dictionary entries
Ntheta    = size(A,3);
% number of snapshots in the data covariance
Nsnapshot = size(Y,3);
% noise power initialization
sigc      = ones(Nfreq,1) * options.noisepower.guess;
% posterior
x_post    = zeros(Nfreq, Ntheta, Nsnapshot);
% minimum (global) gamma
gmin_global   = realmax;
% L1 error
errornorm    = zeros(options.convergence.maxiter,1);

% initialize equal and uncorrelated weights
gamma        = 1*ones(Ntheta,1);
gamma_num    = zeros(Nfreq , Ntheta);
gamma_denum  = zeros(Nfreq , Ntheta);

% Sample Covariance Matrix
SCM = zeros( Nfreq , Nsensor , Nsensor );

for i_f = 1 : Nfreq
    SCM(i_f,:,:) = squeeze2(Y(i_f,:,:)) * squeeze2(Y(i_f,:,:))' / Nsnapshot;
end

%% Main Loop
Af = squeeze(A(1,:,:));
gamma_a = (diag(Af'*squeeze(SCM(1,:,:))*Af));
gamma_b1 = (diag(Af'*squeeze(SCM(1,:,:))*Af,1));
gamma_b2 = gamma_b1;
display(['SBL version ', options.SBL_v ,' initialized.']);
tic
% figure
beta=options.beta;
for j1 = 1 : options.convergence.maxiter
    
    % for error analysis
    gamma_a_Old = gamma_a;
    gamma_b1_Old = gamma_b1;
     gamma_b2_Old = gamma_b2;
    %% gamma update
    
    for i_f = 1 : Nfreq
        
        Af = squeeze(A(i_f,:,:));
        
        % for tridiagonal matrices
        
        %         gamma_a = diag(Gamma);
        %         gamma_b1 = diag(Gamma,-1);
        %         gamma_b2 = diag(Gamma,1);
        
        Afh = Af';
        %SigmaYinv  = inv(sigc(i_f) * eye(Nsensor) + Af * (repmat(gamma_a, [1 Nsensor] ) .* Afh));
        
        SigmaYinv  = inv(sigc(i_f) * eye(Nsensor) + Af * (repmat(gamma_a, [1 Nsensor] ) .* Afh +...
            beta*repmat([gamma_b1;0], [1 Nsensor]) .* [Afh(2:end,:); zeros(1,Nsensor)]+...
            beta*repmat([0;gamma_b2], [1 Nsensor]) .* [zeros(1,Nsensor);Afh(1:end-1,:)]));
%              imagesc(abs(sigc(i_f) * eye(Nsensor) + Af * (repmat(gamma, [1 Nsensor] ) .* Af')));
%         
%                 drawnow
        Yf = squeeze2(Y(i_f,:,:));
        
         ApSigmaYinv  = Af'*SigmaYinv;

        
        % Sum over snapshots and normalize, abs for roundoff errors
        gamma_num(i_f,:)   = sum ( abs ( ( ApSigmaYinv * squeeze2(Y(i_f,:,:)) ).^2 ),2 ) / Nsnapshot;
        
        % positive def quantity, abs for roundoff errors
        gamma_denum(i_f,:) = abs( sum  ( ApSigmaYinv.' .* Af, 1 ) );
        
        U = SigmaYinv*(Yf*Yf')*SigmaYinv;
        V = Nsnapshot*SigmaYinv;
        
        for i = 1: Ntheta-1
            gamma_a(i) = gamma_a_Old(i) * abs(((Af(:,i)'*U*Af(:,i))/(Af(:,i)'*V*Af(:,i)))).^(1/options.fixedpoint);

            gamma_b1(i)= gamma_b1_Old(i)* ((Af(:,i)'*U*Af(:,i+1))/(Af(:,i)'*V*Af(:,i+1))).^(1/options.fixedpoint);
            gamma_b2(i)= gamma_b2_Old(i)* ((Af(:,i+1)'*U*Af(:,i))/(Af(:,i+1)'*V*Af(:,i))).^(1/options.fixedpoint);

        end
        gamma_a(Ntheta) = gamma_a_Old(Ntheta) * (((Af(:,Ntheta)'*U*Af(:,Ntheta))/(Af(:,Ntheta)'*V*Af(:,Ntheta)))).^(1/options.fixedpoint);
%         gamma_b2 = gamma_b1;
        
    end
    gamma = diag(gamma_a) +diag(gamma_b1,1)+diag(gamma_b2,-1);
%     plot(abs(gamma_a))
%     hold on
%     plot(abs(gamma_b1))
% legend('a','b')
% % hold off 
% % % imagesc(abs(gamma))
%     drawnow
%     
%     %% sigma and L2 error using unbiased posterior update
%     aaaa= sum(abs(gamma_a))
%     bbbb=sum(abs(gamma_b1))
    
    %% sigma and L2 error using unbiased posterior update
    
   % locate same peaks for all frequencies
    [ ~ , Ilocs] = findpeaks(abs(gamma_a),'SORTSTR','descend','NPEAKS',Nsource);
    Apeak      = A(:,:,Ilocs);
    
    for i_f = 1 : Nfreq
        
        % only active replicas
        Am     = squeeze2(Apeak(i_f,:,:));
        
        % noise estimate
        sigc(i_f) = real(trace( (eye(Nsensor)-Am*pinv(Am)) * squeeze(SCM(i_f,:,:)) ) / ( Nsensor- Nsource ) );
        
    end

%     
%     plot(gamma_a)
%     title('SBL_Tri')
%     drawnow    


%     for i_f = 1 : Nfreq
%    
%         % noise estimate
%         sigc(i_f) = real(trace(SigmaYinv));
% 
%     end
    
    %% Convergance
    % checks convergance and displays status reports
    
    % convergance indicator
    errornorm(j1) = norm ( gamma_a - gamma_a_Old, 1 ) / norm ( gamma_a, 1 );
    iteration_L1   = j1;
    % look into the past and find best error since then
    if j1 > options.convergence.min_iteration  &&  errornorm(j1) < gmin_global
        gmin_global = errornorm(j1);
        gamma_min   = gamma_a;
        iteration_L1   = j1;
    end
    gamma_min   = gamma_a;
    % inline convergence code
    if j1 > options.convergence.min_iteration && ( errornorm(j1) < options.convergence.error  || iteration_L1 + options.convergence.delay <= j1)
        
        if options.flag == 1
            display(['Solution converged. Iteration: ',num2str(sprintf('%.4u',j1)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])
        end
        break; % goodbye
        
        % not convereged
    elseif j1 == options.convergence.maxiter
        warning(['Solution not converged. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])
        
        % status report
    elseif j1 ~= options.convergence.maxiter  && options.flag == 1 && mod(j1,options.status_report) == 0 % Iteration reporting
        display(['Iteration: ',num2str(sprintf('%.4u',j1)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.' ])
        
    end
    %     gamma_a
    
end

%% Posterior distribution
% x_post - posterior unbiased mean
for i_f = 1 : Nfreq
    
    x_post(i_f,:,:) = (repmat(gamma_a, [1 Nsnapshot] ) +beta*repmat([gamma_b1;0], [1 Nsnapshot]) +beta*repmat([0;gamma_b2], [1 Nsnapshot])) .*...
        (Af' / (sigc(i_f) * eye(Nsensor) +...
        Af * (repmat(gamma_a, [1 Nsensor] ) .* Afh +beta*repmat([gamma_b1;0], [1 Nsensor]) .* [Afh(2:end,:); zeros(1,Nsensor)]+beta*repmat([0;gamma_b2], [1 Nsensor]) .* [zeros(1,Nsensor);Afh(1:end-1,:)])) * squeeze2(Y(i_f,:,:)));
    
    
%     x_post(i_f,:,:) = repmat(gamma_a, [1 Nsnapshot] ) .*...
%         (Af' / (sigc(i_f) * eye(Nsensor) +...
%         Af * (repmat(gamma_a, [1 Nsensor] ) .* Af')) * squeeze2(Y(i_f,:,:)));
    %     x_post(i_f,:,:) = repmat(gamma_a, [1 Nsnapshot] ) .*...
    %         (Af' / (sigc(i_f) * eye(Nsensor) +...
    %         Af * (repmat(gamma_a, [1 Nsensor] ) .* Af')) * squeeze2(Y(i_f,:,:)));
    
    %     imagesc(mean(real(reshape(squeeze(x_post),101,101,size(x_post,3))),3))
    %     drawnow;
    %
    
    gamma = diag(gamma_a)+diag(gamma_b1,-1)+diag(gamma_b1,1);
end

%% function return
% Globla minimum
gamma_a = gamma_min;

% Report section

% vectors containing errors
report.results.error    = errornorm;

% Error when minimum was obtained
report.results.iteration_L1 = iteration_L1;

% General info
report.results.final_iteration.iteration = j1;
report.results.final_iteration.noisepower = sigc;

if options.tic == 1
    report.results.toc = toc;
else
    report.results.toc = 0;
end

% data
report.results.final_iteration.gamma  = gamma_a  ;
report.results.final_iteration.x_post = x_post ;

report.options = options;

end

function b = squeeze2(a)
%   SQUEEZE2
%   Just as squeeze but with transpose to accomodate single snapshot case.
%   This is required because matlab does not allow a singleton dimension
%   at the 'end', e.g., 3x5x1.
%   The entire code might alternatively be re-written to have the number
%   of sensors n~=1 at the end (making this fix obsolete).

if ~ismatrix(a)
    siz = size(a);
    siz(siz==1) = []; % Remove singleton dimensions.
    siz = [siz ones(1,2-length(siz))]; % Make sure siz is at least 2-D
    b = reshape(a,siz);
else
    b = a.';
end

end
