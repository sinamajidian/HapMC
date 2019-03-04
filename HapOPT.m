%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Haplotyping using HapOPT
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


function h=HapOPT(R_used,k)

%k=1; For all heterzygous case, set it to one.
[X, S_opt, Y, ~] = OptSpace(R_used,k,300,.001); %  matrix, rank,number iter, toleranc
X_opt=X*S_opt*Y';
A=X_opt';
[~,colind] = rref(A);
Xsub = A(:, colind(1:k));
h=Xsub'>0;
