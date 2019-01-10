

function [idx_i,idx_j,R_block]=first_block_extractor(R_find)
%%%%%%%%%%%%%%%%%%%%%
% This code seperate the first component of a read matrix which has overlap
% As a part of HapMC

% Input: a matrix that each rows of it is a read {1,-1}
% output: the matrix of first block
%         index of ending element

%Sina Majidian Dec 2018
%Iran University of Science and Technology
%%%%%%%%%%%%%%%%%%%%%



[N_row,N_col]=size(R_find);

if N_row>1 && N_col >1
    for i=2:N_row
        for j=2:N_col
            if sum(abs(R_find(i,1:j-1)))+sum(abs(R_find(1:i-1,j)))==0
                
                bot_left=R_find(i:end,1:j-1);
                up_right=R_find(1:i-1,j:end);
                %For optimizaing purpose, first do it just for a column and a row
                % if sum(abs(R_find(i,1:j-1)))+sum(abs(R_find(1:i-1,j)))==0
                if sum(abs(bot_left(:)))+sum(abs(up_right(:)))==0
                    idx_i=i-1;
                    idx_j=j-1;
                    R_block=R_find(1:idx_i,1:idx_j);
                    return
                end
            end
        end
    end
    idx_i=i;
    idx_j=j;
    R_block=R_find;
else %if N_col==1
    idx_i=[];
    idx_j=[];
    R_block=[];
end
end



