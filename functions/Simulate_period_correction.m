%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function simulate a sequence based on the model proposed
% by Schulze et. al. 2005, an equivalent formula was proposed by 
% Jacoby and Repp 2012 and is used here.
%
% input:
% interstimulus intervals S
%note that S is unconventional since S(t)=OS(t)-OS(t-1) 
% and not like Vorberg and Schulze 2002 where C(t)=OS(t+1)-OS(t)=S(t+1)
% alpha, and beta of the model of Schulze et. al 2005
% sigmaM, sigmaT - motor and internal time-keeper standard deviations
% respectively. 
%
% The model takes the form (proof in appendix II in Jacoby and Repp 2012):
% R(n+1)= (-alpha –beta)A(n) + alpha*A(n-1) + R(n) + T(n)-T(n-1) + M(n)- 2M(n-1)+ M(n-2)
%
% an equivalent formula is:
% S(n+1)+A(n+1) = (2-alpha –beta)A(n) + (alpha-1)A(n-1) + S(n) + T(n)-T(n-1) + M(n)- 2M(n-1)+ M(n-2)
%
% output:
% OSm,ORm - modeled stimulus,repsomse onsets (msec)
%
% Rm-  modeled interesponse interval Rm(t)=ORm(t)-ORm(t-1) 
% note that this is unconventional see remark about s(t) before
% 
% Am - modeled asynchrony: Am(t)=ORm(t)-OSm(t)
% 
% Written by Nori Jacoby, for help please contact: nori.viola@gmail.com
%
%[Schulze et. al. 2012 ] Schulze H-H, Cordes A, Vorberg D (2005) Keeping synchrony while
%tempo changes: accelerando and ritardando. Mus Percept 22:461–477
%
% [Jacoby and Repp 2012]  Nori and Bruno H. Repp.
% “A General Linear Framework for the Comparison and Evaluation of Models of Sensorimotor Synchronization.” 
% Biological Cybernetics, vol. 106/3, 135-154, DOI: 10.1007/s00422-012-0482-x.
%
%   For more information see: Nori Jacoby, Naftali Tishby, Bruno H. Repp, Merav Ahissar and Peter E. Keller (2015)
%
%
% =====  CODE BY: Nori Jacoby (nori.viola@gmail.com)
%
% NOTE!!!!!  this simulation outputs Tm the timekeeper!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==========================================================================
% ===== For information please contact: Nori Jacoby 
% ===== Nori.viola@gmail.com
% =====
% ===== If you are using the code,
% ===== Please cite this version:
% ===== 
% ===== Jacoby, Nori, Peter Keller, Bruno H. Repp, Merav Ahissar and Naftali Tishby. 
% ===== "Parameter Estimation of Linear Sensorimotor Synchronization Models: Phase Correction, Period Correction and Ensemble Synchronization." 
% ===== Special issue of Timing & Time Perception (RPPW).
% ==========================================================================

function  [Rm,Am,OSm,ORm,Tm] = Simulate_period_correction(S,alpha,beta,sigmaT,sigmaM)

assert(sum(S<=0)==0);

N=length(S);
OSm=cumsum(S);

Am=zeros(N,1);
Rm=zeros(N,1);
ORm=zeros(N,1);


Tn=sigmaT*randn(N+1,1);
Mn=sigmaM*randn(N+2,1);
Hn=Tn(2:end)-Tn(1:(end-1)) + Mn(3:end)- 2*Mn(2:(end-1)) + Mn(1:(end-2));

ORm(1:2)=OSm(1:2);
Rm(1:2)=S(1:2);
Am(1:2)=OSm(1:2)-ORm(1:2);
Tm(1:2)=S(1:2);

for K=(2):(N-1),
    ORm(K+1)=ORm(K)+Hn(K) + (-alpha-beta)*Am(K)+ (alpha)*Am(K-1) + Rm(K);
    Am(K+1)=ORm(K+1)-OSm(K+1);
    Rm(K+1)=ORm(K+1)- ORm(K);
    Tm(K+1)=Am(K+1)+S(K+1)-(1-alpha)*Am(K)-Hn(K);
end

