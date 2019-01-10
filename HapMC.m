%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Haplotyping using HapSVT

% Input: read matrix in .mat format
% output: a text file containg haplotypes and corresponding variant index

% the longest block is considered 2000
%Sina Majidian Dec 2018
%Iran University of Science and Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clearvars
% close all
load('data_geraci/mat_used/e01c10.mat') %/SinaMc/University/codeMy/CodePaper1-UsingMatrixcompletion/Senario_SimulatedData/data_geraci/mat_used
data_number=100;



rr_o=zeros(1,data_number);
rr_a=zeros(1,data_number);

swe_o=zeros(1,data_number);
swe_a=zeros(1,data_number);


NumberRowsofMatrix=NumberRowsofMatrix(:,2);
NumRow=[0;cumsum(NumberRowsofMatrix)];

for jj=6:data_number




% rr_svt=zeros(1,Data_number);rr_alt=zeros(1,Data_number); rr_svt=zeros(1,Data_number);
% swe_svt=zeros(1,Data_number);swe_alt=zeros(1,Data_number);swe_svt=zeros(1,Data_number);
% t_svt=zeros(1,Data_number);t_alt=zeros(1,Data_number); t_svt=zeros(1,Data_number);
% MEC_svt=zeros(1,Data_number);MEC_alt=zeros(1,Data_number); MEC_svt=zeros(1,Data_number);

% for jj=1:Data_number
%  for jj=1:Data_number


R_all=ReadMatrixe0(NumRow(jj)+1:NumRow(jj+1),:);
H_ex_all=Haplotypese0([2*jj,2*jj-1],:);

% R_all=ReadMatrix(NumRow(jj)+1:NumRow(jj+1),:);
% H_ex_all=Haplotypes1([2*jj,2*jj-1],:);


%%% the data is not all heter  !
indx_het=[];
for j=1:size(H_ex_all,2)
    if H_ex_all(1,j)~=H_ex_all(2,j)
    indx_het=[indx_het, j];
    end
end
R_used_col=R_all(:,indx_het);
H_ex=H_ex_all(:,indx_het);
nonzeor_idx_row=find(sum(abs(R_used_col'))>1); % those rows with at least two nonzero
R_used=R_used_col(nonzeor_idx_row,:);

% R_used=R_all(1:59,1:30);
% h_exct1a=h_exct1(1:30);
% sum(abs(R))



%%%% haplotpying  using OPT

% [X,S_opt,Y,dist] = OptSpace(R_used,1,500,.00001); %(R,[],500,.001)  matrix, rank,number iter, toleranc

% X_opt=X*S_opt*Y';
% A=X_opt';

%         %%%%% haplotpying  using SVT
        omega=find(R_used);
        [U,S,V,numiter]= FPC(size(R_used),omega,R_used(omega),1);  %FPC(n,Omega,b,mu_final,maxiter,tol)
        X1=U*S*V';[u_sv,S_sv,v_sv]=svds(X1,1); X_svt=u_sv*S_sv*v_sv';
        A=X_svt';

%         %%%%% haplotpying  using Nuclear minimization
%         full(R_used)
% %         omega=find(R_used);
%         [N,l] = size(R_used);
% 
%         cvx_begin quiet
%         variable X(N,l)
%         minimize norm_nuc(X)
%         subject to
%         norm(X(omega)-R_used(omega),2)<=.1;
%         cvx_end
% 
%         X1=full(X);  % converts from sparse
%         [u,S,v]=svds(X1,1); X_nc=u*S*v';
%         A=X_nc';


[~,colind] = rref(A);
Xsub = A(:, colind(1:1));
h=2*(Xsub'>0)-1;
H_o=[h;-h];
[rr_sv,swe_sv,sh_sv,l_sv] = statcal(H_ex,H_o);

% [rr_nc,swe_nc,sh_nc,l_nc]

%%%% althap
% tic
% H_alt= AltHapedit(R_used);
% toc
% [rr_alt,swe_alt,sh_alt,l_alt] = statcal(H_ex,H_alt);



rr_o(jj)=rr_sv
% rr_a(jj)=rr_alt

swe_o(jj)=swe_sv;
% swe_a(jj)=swe_alt;


% [rr_svt,swe_svt,sh_svt,l_svt]
% [rr_alt,swe_alt,sh_alt,l_alt]

% h_exct1a=h_exct1;
% sum(h_alt+h_exct1a')/2
% min(abs([sum(h(2,:)+h_exct1a),sum(h(2,:)-h_exct1a),sum(h(1,:)-h_exct1a),sum(h(1,:)+h_exct1a)]))/1
% 
% 
% min(abs(sum(h(1,:)+h1)),abs(sum(h(1,:)-h1)))
% 

end


save('out/e02c10_out_nc.mat','-v7.3')
