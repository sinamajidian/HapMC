%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Haplotyping using HapSVT
%
% Input: Read matrix in which all reads have overlaps.
% output: A matrix or vector of haplotype (witout positions)
%
%CVX should be installed beforehand. (http://cvxr.com/cvx/)
%
% This code is part of HapMC package.
%
%Sina Majidian Dec 2018
%Iran University of Science and Technology
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=HapNuc(R_used,k)

omega=find(R_used);
[N, l] = size(R_used);

cvx_begin quiet
variable X(N,l)
minimize norm_nuc(X)
subject to
norm(X(omega)-R_used(omega),2)<=.1;
cvx_end

X1=full(X);  % converts from sparse
[u,S,v]=svds(X1,k); X_nc=u*S*v';
A=X_nc';
[~,colind] = rref(A);
Xsub = A(:, colind(1:k));
h=Xsub'>0;