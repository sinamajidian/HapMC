%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Haplotyping using HapSVT
%
% Input: Read matrix in which all reads have overlaps.
% output: A matrix or vector of haplotype (witout positions)
%
% This code is part of HapMC package.
%
%
%Sina Majidian Dec 2018
%Iran University of Science and Technology
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function h=HapSVT(R_used,k)

% one block of overlapped reads




%%%%%%%% haplotpying  using SVT
omega=find(R_used);
[U, S, V, numiter]= SVT(size(R_used),omega,R_used(omega),.001);%FPC(n,Omega,b,mu_final,maxiter,tol)
X1=U*S*V';
[u_sv,S_sv,v_sv]=svds(X1,1);

X_svt=u_sv*S_sv*v_sv';
[~,colind] = rref(X_svt);
Xsub = X_svt(:, colind(1:k));
h=Xsub'>0;


