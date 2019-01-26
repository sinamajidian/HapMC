function [error]=HapMC(R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Haplotyping using HapSVT, HapOPT, HapNuc

% Input: read matrix in .mat format
% output: a text file containg haplotypes and corresponding variant index

%?he longest block is considered 2000
%Sina Majidian Dec 2018
%Iran University of Science and Technology


% OptSpace.m and the SVT.m should be downloaded separately.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_in=R;

fileID_opt = fopen('chr1_opt.hap','w'); % The output file
fileID_nuc = fopen('chr1_opt.hap','w'); % The output file
fileID_svt = fopen('chr1_opt.hap','w'); % The output file



while columnNumber_blocks_accm(end)<size(R_in,2)-1  % run for each block until last column
    block_idx = block_idx+1;
    if size(R_in,2)<columnNumber_blocks_accm(end)+2000 | size(R_in,1)<rowNumber_blocks_acc(end)+2000
        % the start column is the last value of accumaltion all haploblock length +1
        % the end column is .+2000
        R_sliced2000=R_in(rowNumber_blocks_acc(end)+1:end,columnNumber_blocks_accm(end)+1:end);
    else
        R_sliced2000=R_in(rowNumber_blocks_acc(end)+1:rowNumber_blocks_acc(end)+2000,columnNumber_blocks_accm(end)+1:columnNumber_blocks_accm(end)+2000);
    end
    if start==1  % for the first block
        start=0;rowNumber_blocks_acc=[];columnNumber_blocks_accm=[];
    end
    %Exctracting each block of read matrix which all reads have overlap
    [rowNumber_block,columnNumber_block,R_block]=first_block_extractor(R_sliced2000);%block_length is the haploBlock length
    rowNumber_blocks=[rowNumber_blocks, rowNumber_block];columnNumber_blocks=[columnNumber_blocks, columnNumber_block];
    rowNumber_blocks_acc=cumsum(rowNumber_blocks);columnNumber_blocks_accm=cumsum(columnNumber_blocks);
    if  size(R_block,2)>1
        %processing of read matrix block
        nonzeor_idx_row=find(sum(abs(R_block'))>1); % those rows with at least two nonzero
        R_used=R_block(nonzeor_idx_row,:);
        
        %%%% haplotpying  using OPT
        k=2; % For all het case, set it to one.
        [X,S_opt,Y,dist] = OptSpace(R_used,k,300,.001); %  matrix, rank,number iter, toleranc
        X_opt=X*S_opt*Y';
        A=X_opt';
        [~,colind] = rref(A);
        Xsub = A(:, colind(1:1));
        h=Xsub'>0;
        if block_idx==1 %for reporting the index of estimated allele of variant to not lossing
            %the index in overall, we remove some columns
            %columnNumber_blocks_accm:  these number are after removing empty columns
            %indces_block:  these number is before removing empty columns
            indces_block=nonzero_col_idx(1:columnNumber_blocks_accm(1));
        else
            indces_block=nonzero_col_idx(columnNumber_blocks_accm(block_idx-1)+1:columnNumber_blocks_accm(block_idx));
        end
        H_with_ind=[indces_block',h'];
        fprintf(fileID_opt,'BLOCK\t%d\t%d\n',indces_block(1),indces_block(end));
        fprintf(fileID_opt,'%d\t%d\n',H_with_ind');
        
        
        
        %%%%%%%% haplotpying  using SVT
        omega=find(R_used);
        [U,S,V,numiter]= SVT(size(R_used),omega,R_used(omega),.001);  %FPC(n,Omega,b,mu_final,maxiter,tol)
        X1=U*S*V';[u_sv,S_sv,v_sv]=svds(X1,1); X_svt=u_sv*S_sv*v_sv';
        A=X_svt';
        A=X_opt';
        [~,colind] = rref(A);
        Xsub = A(:, colind(1:k));
        h=Xsub'>0;
        if block_idx==1 %for reporting the index of estimated allele of variant to not lossing
            indces_block=nonzero_col_idx(1:columnNumber_blocks_accm(1));
        else
            indces_block=nonzero_col_idx(columnNumber_blocks_accm(block_idx-1)+1:columnNumber_blocks_accm(block_idx));
        end
        H_with_ind=[indces_block',h'];
        fprintf(fileID_svt,'BLOCK\t%d\t%d\n',indces_block(1),indces_block(end));
        fprintf(fileID_svt,'%d\t%d\n',H_with_ind');
        
        
        
        %%%%% haplotpying  using Nuclear minimization
        full(R_used)
        omega=find(R_used);
        [N,l] = size(R_used);
        
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
        if block_idx==1 %for reporting the index of estimated allele of variant to not lossing
            indces_block=nonzero_col_idx(1:columnNumber_blocks_accm(1));
        else
            indces_block=nonzero_col_idx(columnNumber_blocks_accm(block_idx-1)+1:columnNumber_blocks_accm(block_idx));
        end
        H_with_ind=[indces_block',h'];
        fprintf(fileID_nuc,'BLOCK\t%d\t%d\n',indces_block(1),indces_block(end));
        fprintf(fileID_nuc,'%d\t%d\n',H_with_ind');
        
    end
end
fclose(fileID);


save('data.mat','-v7.3')

error = 0;


%function [X S Y dist] = OptSpace(M_E,r,niter,tol)
% An algorithm for Matrix Reconstruction from a partially revealed set. 
% See "Matrix Completion from a Few Entries"(http://arxiv.org/pdf/0901.3150) for details
% Usage :
% [X S Y dist] = OptSpace(A,r,niter,tol);
% [X S Y dist] = OptSpace(A);
% 
% INPUT :
% A     :  The partially revealed matrix.
%          Sparse matrix with zeroes at the unrevealed indices.
%
% r     :  The rank to be used for reconstruction. Use [] to guess the rank.
% niter :  The max. no. of iterations. Use [] to use default (50).
% tol   :  Stop iterations if norm( (XSY' - M_E).*E , 'fro' )/sqrt(|E|) < tol, where
%        - E_{ij} = 1 if M_{ij} is revealed and zero otherwise, 
%        - |E| is the size of the revealed set.				
%        - Use [] to use the default (1e-6)
%
%
% OUTPUT :
% X      : A size(A,1)xr matrix
% S      : An rxr matrix
% Y      : A size(A,2)xr matrix
% such that M_hat = X*S*Y' 
% dist   : A vector containing norm( (XSY' - M_E).*E , 'fro' )/sqrt(|E|) at each
%          successive iteration
%
% Date : 21st April, 2009
% COPYRIGHT 2009 Raghunandan H. Keshavan, Andrea Montanari, Sewoong Oh





%function [U,Sigma,V,numiter,out]  = SVT(n,Omega,b,tau,delta,maxiter,tol,EPS)
% [U,Sigma,V,numiter,output]  = SVT(n,Omega,b,tau,delta,maxiter,tol,EPS)


%http://svt.stanford.edu/code/

% Finds the minimum of   tau ||X||_* + .5 || X ||_F^2 
%
% subject to P_Omega(X) = P_Omega(M)
%
% using linear Bregman iterations
%
% Usage:  [U,S,V,numiter]  = SVT(n,Omega,b,delta,maxiter,tol)
%
% Inputs:
%
%   n - size of the matrix X assumed n(1) by n(2). If n is a single integer, it
% is understood that n(1) = n(2). 
%
%   Omega - set of observed entries.  Should be linearly indexed.
%
%   b - data vector of the form M(Omega)
%
%   tau - parameter defining the objective functional 
%
%   delta - step size.  Choose delta less than 2 to be safe but
%       conservative; choose delta closer to n(1)*n(2)/length(Omega)
%       to be riskier (i.e. algorithm may diverge)
%
%   maxiter - maximum number of iterations
%
%   tol - stopping criteria (default: 1e-4)
%
%   EPS - noise constraint.  This relaxes the constraints, so that they
%       are now of the form | X(i,j) - M(i,j) | <= EPS,
%       for all indices (i,j) in omega.  Default: 0
%
% Outputs: matrix X stored in SVD format X = U*diag(S)*V'
% 
%   U - n1xr left singular vectors 
% 
%   S - rx1 singular values
%
%   V - n2xr right singular vectors 
%
%   numiter - number of iterations to achieve convergence
%
%   output - a structure with data from each iteration.  Includes:
%       output.nuclearNorm  - nuclear norm of current iterate
%       output.rank         - rank of current iterate
%       output.time         - time taken for one iteraration
%       output.residual     - the relative residual, norm(x-b)/norm(b)
% Description: 
% Reference:
%
%    Cai, Candes and Shen
%    A singular value thresholding algorithm for matrix completion
%    Submitted for publication, October 2008.
%
%    See also more general code as part of the TFOCS package,
%    available at tfocs.stanford.edu as of November 2010.
%
% Written by: Emmanuel Candes
% Email: emmanuel@acm.caltech.edu
% Created: October 2008
% Efficient mex-file and PROPACK version: Stephen Becker, Nov 2008
% Modified: Stephen Becker, March 2009
% Modified: Stephen Becker, May 2009  works with complex numbers
% Modified: Farshad Harirchi and Stephen Becker, April 2011, fixing a bug.

%function [U,S,V,numiter]  = FPC(n,Omega,b,mu_final,maxiter,tol)
% [U,Sigma,V,numiter]  = FPC(n,Omega,b,mu_final,maxiter,gtol)
%
% Finds mininum  mu ||X||_* + 1/2 || A(X) - b ||_2^2
%
%   where A(X) is the projection of X onto the set of indices
%   in Omega.
%
% For efficiency, the algorithm uses continuation (i.e. a series of
%   mu, aka the "outer loop"), until mu = mu_final.
%
% maxiter controls maximum number of iterations per inner loop
%
% Outputs:
%   U,V and Sigma are singular vectors and singular values, where
%       X = U*Sigma*V'
%   numiter is the number of iterations over all inner and outer loops

% Reference:
%
%   "Fixed point and Bregman iterative methods for matrix
%       rank minimization."
%   Shiquian Ma, Donald Goldfarb, Lifeng Chen, October 2008
%   ftp://ftp.math.ucla.edu/pub/camreport/cam08-78.pdf

% code by Stephen Becker, srbecker@caltech.edu, March 2009

% May 2009: adding support for complex matrices


