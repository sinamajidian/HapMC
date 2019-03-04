%function [error]=HapMC(fragment_file,Hap_algorithm)
fragment_file='data/fragment_sample.txt';
Hap_algorithm='O';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Haplotyping using HapSVT, HapOPT, HapNuc
%
% Input: read matrix in .mat format
% algorithm canbe either 'HapOPT' or 'HapSVT' or 'HapNuc'.
% output: a text file containg haplotypes and corresponding variant index
%?he longest block is considered 2000
%
%
%Sina Majidian Dec 2018
%Iran University of Science and Technology
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


R_in=convert_frag_mat(fragment_file,SNP_Numbers);


fileID_opt = fopen('Reconstructed_Haplotype.txt','w'); % The output file



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
        
        if Hap_algorithm=='O'
            
            h=HapOPT(R_used,k);
            
        elseif Hap_algorithm=='S'
            h=HapSVT(R_used,k);
            
        elseif Hap_algorithm=='N'
            h=HapNuc(R_used,k);
            
            
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
    
    
    
