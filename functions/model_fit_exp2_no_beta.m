% bGLS without period correction based on the difference equation, to
% compare with the period correction model fit
% Missing \ excluded taps are nan in As and R (both the interval after and the interval before)
% then we simply erase the relevant rows and columns in the covariance matrix of the noise - 
% so all relevant connections between variables are just part of the unmodelled noise.
%
% For more information about the method see:
% Jacoby, Tishby, Repp, Ahissar and Keller (2015)
% and Vishne, Jacoby et al. (2021)
% Please cite if you use this code
%
% Let OR(t) be the response onset at time t
% Let OS(t) be the stimulus onset at time t
% Let R(t) be the interesponse interval R(t)=OR(t)-OR(t-1)
% Let A(t) be the asynchornies A(t)=OR(t)-OS(t)
% note that R is is slightly diffrent than the notation of 
% Vorberg and Shultze 2002 where:
% I(t)=OR(t+1)-OR(t)=R(t+1)
% the empirical means should be computed outside.
% 
% CODE BY: Nori Jacoby (nori.viola@gmail.com) & Gal Vishne (gal.vishne@gmail.com)

function [alphas, st, sm, logp]=model_fit_exp2_no_beta(R,As,MEAN_A)

% parameters of GLS method
ITER=20; % number of iterations
TRESH=1e-3; %when changes betwen iterations are smaller than this, stop!

N=size(R,1)-2; % number of samples
P=size(As,2);  % number of partners
assert(size(R,1)==size(As,1));

for p=1:P
     As(:,p)=As(:,p)-MEAN_A(p);
end

b3=R(3:end)-R(2:end-1);         % create the matrices
A3=[As(2:end-1,:)-As(1:end-2,:)];
    
missing_inds = any(isnan([b3,A3]),2);
b3(missing_inds)=[];
A3(missing_inds,:)=[];

% init acvf
K11=2;
K12=-1;
K13=0;

zold=zeros(2*P,1)-9999; %init to invalid value

% do the BGLS iterations
for iter=1:ITER
        CC=diag(K11*ones(1,N),0)+ diag(K12*ones(1,N-1),1)+ diag(K12*ones(1,N-1),-1)+...
            diag(K13*ones(1,N-2),2)+ diag(K13*ones(1,N-2),-2);
        CC(missing_inds,:)=[];CC(:,missing_inds)=[];
            z=((A3')*(CC\A3))\((A3')*(CC\b3)); % compute GLS
        d=A3*z-b3; % compute residual noise
        
        K=cov(d(1:(end-1)),d(2:end));    %estimate residual acvf
        K11=(K(1,1)+K(2,2))/2;
        K12=K(1,2);

        % apply bound
        if K12 < (-5/8)*K11
            K12=(-5/8)*K11;
        end
        if K12 > (-1/2)*K11
            K12 = (-1/2)*K11;
        end
        K13 = -0.5*(K11+2*K12);
        % if allready obtain local maxima there is no point to continue...
        if (sum(abs(z-zold))<TRESH)
            break;
        end
        zold=z;
end % end of bGLS iterations.
    
%output the values
alphas=-z(1:P);
sm=sqrt(-0.5*(K11+2*K12));
st=sqrt(2*K11+3*K12);
    
% log likelihood
d=d';
const = -0.5 * length(d) * log(2*pi);
term1 = -0.5 * sum((d / CC) .* d, 2); % N x 1
term2 = const - sum(log(diag(chol(CC))));    % scalar
logp = term1' + term2;