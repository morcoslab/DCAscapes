function [c_average,h_average] = masterDCAparameters_length(inputfile,len)

%%%% unfixed length 
% Direct Coupling Analysis (DCA) Scapes - parameters estimation
% couplings and local fields parameters estimation
%
% INPUTS:
%   inputfile  - file containing the FASTA alignment
%	len 	   -  the length of the sequence
%
%
% OUTPUTS:
%   h_average 	    - q x N matrix with the average local fields obtained
%                       from the ones calculated pairwise.
%   c_average       - (q*N) x (q*N) matrix with the eij(alpha,beta)
%                     couplings. Each submatrix i,j (size q x q) contains
%                     the couplings for the i,j pair.
%
%
%
%
% This software and accompanying documents are implementated based on the 2011 DCA paper
% (F Morcos, A Pagnani & B Lunt et al, 2011 ) and code by :
%             2011/12 - Andrea Pagnani and Martin Weigt
%                       andrea.pagnani@gmail.com
%                       martin.weigt@upmc.fr
%
% This implementation and accompanying scripts (DCAparameters*.m and Fastahamiltonian.m, newdca.m)
% include changes to process SEQRS (RNA) data and calculation of Hamiltonians and
% other metrics to study Protein-RNA interactions
%
% Copyright for this implementation:
%
%             2018/4  - Qin Zhou, José Alberto De la Paz and Faruck Morcos
%                        qxz142130@utdallas.edu
%                        alberto_nay@hotmail.com
%                        faruckm@utdallas.edu
%
%
% Any publication resulting from applications of DCA and DCA-scapes should cite:
%
%
%     Q Zhou, N Kunder, José Alberto De la Paz, AE. Lasley, VD.Bhat,
%     F Morcos, ZT. Campbell (2018),Global pairwise RNA interaction
%     landscapes reveal corefeatures of protein recognition.
%
%     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
%     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
%     analysis of residue co-evolution captures native contacts across
%     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
%
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=cell(1,4);
familycouplings=cell(1,4);
align=cell(1,4);
for i=1:4
	[h{i},familycouplings{i},~]=DCAparameters_1(inputfile,2,i);
end

[c_average,h_average] = average_couplings_localfields(familycouplings{1},familycouplings{2},familycouplings{3},familycouplings{4},h{1},h{2},h{3},h{4},len);

end


function [c_average,h_average] = average_couplings_localfields(familycouplings1,familycouplings2,familycouplings3,familycouplings4,h1,h2,h3,h4,len)
% Direct Coupling Analysis (DCA) scape - averages the couplings and local fileds parameters from DCAparameters*.m
% function average_couplings_localfields(familycouplings1,familycouplings2,familycouplings3,familycouplings4,h1,h2,h3,h4)
%
% INPUTS:
%   familycouplings1 - output Familycouplings from DCAparameters_1
%   familycouplings2 - output Familycouplings from DCAparameters_2
%   familycouplings3 - output Familycouplings from DCAparameters_3
%   familycouplings4 - output Familycouplings from DCAparameters_4
%
%   h1 - output local fileds (h) from DCAparameters_1
%   h2 - output local fileds (h) from DCAparameters_2
%   h3 - output local fileds (h) from DCAparameters_3
%   h4 - output local fileds (h) from DCAparameters_4
%
% OUTPUTS:
%   c_average	- mean of the couplings for each of the four nucleotide gauged states.
%   h_average	- mean of the local fields for each of the four nucleotide gauged states.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% transfer clouplings 1 to 4
c_4_1=zeros(len*4);
for i=[1:4:len*4]
    for (j=[1:4:len*4])
    c_4_1(i,j)=familycouplings1(i+3,j+3);
    c_4_1(i,j+1)=familycouplings1(i+3,j+1);
    c_4_1(i,j+2)=familycouplings1(i+3,j+2);
    c_4_1(i,j+3)=familycouplings1(i+3,j);

    c_4_1(i+1,j)=familycouplings1(i+1,j+3);
    c_4_1(i+1,j+1)=familycouplings1(i+1,j+1);
    c_4_1(i+1,j+2)=familycouplings1(i+1,j+2);
    c_4_1(i+1,j+3)=familycouplings1(i+1,j);

    c_4_1(i+2,j)=familycouplings1(i+2,j+3);
    c_4_1(i+2,j+1)=familycouplings1(i+2,j+1);
    c_4_1(i+2,j+2)=familycouplings1(i+2,j+2);
    c_4_1(i+2,j+3)=familycouplings1(i+2,j);

    c_4_1(i+3,j)=familycouplings1(i,j+3);
    c_4_1(i+3,j+1)=familycouplings1(i,j+1);
    c_4_1(i+3,j+2)=familycouplings1(i,j+2);
    c_4_1(i+3,j+3)=familycouplings1(i,j);
    end
end


%%% transfer couplings 2 to 4
c_4_2=zeros(len*4);
for i=[1:4:len*4]
    for (j=[1:4:len*4])
    c_4_2(i,j)=familycouplings2(i,j);
    c_4_2(i,j+1)=familycouplings2(i,j+3);
    c_4_2(i,j+2)=familycouplings2(i,j+2);
    c_4_2(i,j+3)=familycouplings2(i,j+1);

    c_4_2(i+1,j)=familycouplings2(i+3,j);
    c_4_2(i+1,j+1)=familycouplings2(i+3,j+3);
    c_4_2(i+1,j+2)=familycouplings2(i+3,j+2);
    c_4_2(i+1,j+3)=familycouplings2(i+3,j+1);

    c_4_2(i+2,j)=familycouplings2(i+2,j);
    c_4_2(i+2,j+1)=familycouplings2(i+2,j+3);
    c_4_2(i+2,j+2)=familycouplings2(i+2,j+2);
    c_4_2(i+2,j+3)=familycouplings2(i+2,j+1);

    c_4_2(i+3,j)=familycouplings2(i+1,j);
    c_4_2(i+3,j+1)=familycouplings2(i+1,j+3);
    c_4_2(i+3,j+2)=familycouplings2(i+1,j+2);
    c_4_2(i+3,j+3)=familycouplings2(i+1,j+1);

    end
end



%%% transfer couplings 2 to 4
c_4_3=zeros(len*4);
for i=[1:4:len*4]
    for (j=[1:4:len*4])
    c_4_3(i,j)=familycouplings3(i,j);
    c_4_3(i,j+1)=familycouplings3(i,j+1);
    c_4_3(i,j+2)=familycouplings3(i,j+3);
    c_4_3(i,j+3)=familycouplings3(i,j+2);

    c_4_3(i+1,j)=familycouplings3(i+1,j);
    c_4_3(i+1,j+1)=familycouplings3(i+1,j+1);
    c_4_3(i+1,j+2)=familycouplings3(i+1,j+3);
    c_4_3(i+1,j+3)=familycouplings3(i+1,j+2);

    c_4_3(i+2,j)=familycouplings3(i+3,j);
    c_4_3(i+2,j+1)=familycouplings3(i+3,j+1);
    c_4_3(i+2,j+2)=familycouplings3(i+3,j+3);
    c_4_3(i+2,j+3)=familycouplings3(i+3,j+2);

    c_4_3(i+3,j)=familycouplings3(i+2,j);
    c_4_3(i+3,j+1)=familycouplings3(i+2,j+1);
    c_4_3(i+3,j+2)=familycouplings3(i+2,j+3);
    c_4_3(i+3,j+3)=familycouplings3(i+2,j+2);
    end
end


c_average=(c_4_1+c_4_2+c_4_3+familycouplings4)/4;



h_4_1=zeros(4,len);
%%% transfer 1 to 4
for i=1:len
    h_4_1(1,i)=h1(4,i);
    h_4_1(2,i)=h1(2,i);
    h_4_1(3,i)=h1(3,i);
    h_4_1(4,i)=h1(1,i);
end

h_4_2=zeros(4,len);
%%% transfer 2 to 4
for i=1:len
    h_4_2(1,i)=h2(1,i);
    h_4_2(2,i)=h2(4,i);
    h_4_2(3,i)=h2(3,i);
    h_4_2(4,i)=h2(1,i);
end

h_4_3=zeros(4,len);
%%% transfer 3 to 4
for i=1:len
    h_4_3(1,i)=h3(1,i);
    h_4_3(2,i)=h3(2,i);
    h_4_3(3,i)=h3(4,i);
    h_4_3(4,i)=h3(3,i);
end

h_average=(h_4_1+h_4_2+h_4_3+h4)/4;

end

function [h,familycouplings,align]=DCAparameters_1(inputfile,stype,rtype)
% Direct Coupling Analysis (DCA) Scapes - parameters estimation
% couplings and local fields parameters estimation, where all the couplings and fields measured relative to nucleotide A were set to 0
%
% function DCAparameters_1(inputfile,stype)
%
% INPUTS:
%   inputfile  - file containing the FASTA alignment
%   stype      - species type:  1 for proteins
%			        2 for RNA and DNA
%   rtype      - relative type: 1 for couplings and fields measured relative to nucleotide A were set to 0 
%				2 for couplings and fields measured relative to nucleotide C were set to 0
%				3 for couplings and fields measured relative to nucleotide G were set to 0
%				4 for couplings and fields measured relative to nucleotide T/U were set to 0
%
% OUTPUTS:
%   h 		    - q x N matrix with the average local fields obtained
%                     from the ones calculated pairwise.
%   Familycouplings - (q*N) x (q*N) matrix with the eij(alpha,beta)
%                     couplings. Each submatrix i,j (size q x q) contains
%                     the couplings for the i,j pair.
%   align - Translated alignment Aminoacid --> Numbers
%
% SOME RELEVANT VARIABLES:
%   N        number of residues in each sequence (no insert)
%   M        number of sequences in the alignment
%   Meff     effective number of sequences after reweighting
%   q        equal to 21 (20 aminoacids + 1 gap) for stype 1
%            equal to 4 for stype 2 (5 if given an alignment with gaps)
%   align    M x N matrix containing the alignmnent
%   Pij_true N x N x q x q matrix containing the reweigthed frequency
%            counts.
%   Pij      N x N x q x q matrix containing the reweighted frequency
%            counts with pseudo counts.
%   C        N(q-1) x N(q-1) matrix Pij_true,Pi_truecontaining the covariance matrix.
%
%   NOTE: The gauge chosen for this implementation is such that local fields and
%	  couplings of q are fixed to be zero. Zeros are explicitly included on the
%	  corresponding matrices.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic
    pseudocount_weight = 0.5; % relative weight of pseudo count
    theta = 0.0;              % threshold for sequence id in reweighting


    [N,M,q,align] = return_alignment(inputfile,stype,rtype);


	if (M>80000)
		 tic
		[Pij_true,Pi_true, Meff]=newCompute_True_Frequencies(align,M,N,q, theta);
		t_stop =toc;
		fprintf ( 1, 'Elapsed CPU time frequencies= %f\n', t_stop );
	else
		tic
		[Pij_true,Pi_true, Meff]=oldCompute_True_Frequencies(align,M,N,q, theta);
		t_stop =toc;
		fprintf ( 1, 'Elapsed CPU time frequencies= %f\n', t_stop );
	end

    fprintf('### N = %d M = %d Meff = %.2f q = %d\n', N,M,Meff,q);
    [Pij,Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q);
    C = Compute_C(Pij,Pi,N,q);
    invC = inv(C);
    familycouplings=nicematrix(-invC,q);
    Pairwisehfield=Compute_Results(Pi, invC, N, q);
    Pairwisehfield=symmetriclocal(Pairwisehfield,N,q);
    [h,~]=averagehfield(Pairwisehfield);
    t_stop= toc;
    fprintf ( 1, 'Elapsed CPU time= %f\n', t_stop);
end

function [N,M,q,Z] = return_alignment(inputfile,stype,rtype)
% reads alignment from inputfile, removes inserts and converts into numbers

    align_full = fastaread(inputfile);
    M = length(align_full);
    ind = align_full(1).Sequence ~= '.' & ...
        align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Z = zeros(M,N);
            for i=1:M
                counter = 0;
                for j=1:length(ind)
                    if( ind(j) )
                        counter = counter + 1;
                        Z(i,counter)=letter2number( align_full(i).Sequence(j), stype, rtype);
                    end
                end
            end
    q = max(max(Z));
end

function Pairwisehfield=Compute_Results(Pi,invC, N,q)
% computes and prints the mutual and direct informations

    Pairwisehfield=zeros(N*q,2*N);
    for i=1:N
        for j=(i+1):N
            % direct information from mean-field
            W_mf = ReturnW(invC,i,j,q);
            Pairwisehfield(((i-1)*q+1):i*q,(2*j-1):2*j)= bp_link(i,j,W_mf,Pi,q);
         end
    end
end

function [Pij_true,Pi_true,Meff] = newCompute_True_Frequencies(align,M,N,q,theta)
% computes reweighted frequency counts

    W = ones(1,M);
    if( theta > 0.0 )
        parfor seq=1:M
            for comparing=1:M
                W(seq)=W(seq)+(pdist([align(seq);align(comparing)],'hamm')<theta);
            end
        end
        W = (1./W);
    end
    Meff=sum(W);

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end

function [Pij_true,Pi_true,Meff] = oldCompute_True_Frequencies(align,M,N,q,theta)
% computes reweighted frequency counts

    W = ones(1,M);
    if( theta > 0.0 )
        W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
    end
    Meff=sum(W);

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end

function x=letter2number(a,stype,rtype)
	switch (stype)
	case 1
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'
            x=2;
        case 'C'
            x=3;
        case 'D'
            x=4;
        case 'E'
            x=5;
        case 'F'
            x=6;
        case 'G'
            x=7;
        case 'H'
            x=8;
        case 'I'
            x=9;
        case 'K'
            x=10;
        case 'L'
            x=11;
        case 'M'
            x=12;
        case 'N'
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
	case 2

	switch(rtype)
	
	case 1
    		switch(a)
	        % full AA alphabet
        	case 'A'
        	     x=4;
        	case 'C'
        	    x=2;
        	case 'G'
        	    x=3;
        	case 'T'
        	    x=1;
        	case 'U'
        	    x=1;
        	case '-'
        	    x=5;
        	otherwise
        	    x=4;  % NOTE:  Correct to x=5;
    		end

	case 2
		switch(a)
        	% full AA alphabet
	        case 'A'
        	     x=1;
        	case 'C'
        	    x=4;
        	case 'G'
        	    x=3;
        	case 'T'
        	    x=2;
        	case 'U'
        	    x=2;
        	case '-'
        	    x=5;
        	otherwise
        	    x=1;  % NOTE:  Correct to x=5;
		end

	case 3
		 switch(a)
	        % full AA alphabet
        	case 'A'
        	     x=1;
        	case 'C'
        	    x=2;
        	case 'G'
        	    x=4;
        	case 'T'
        	    x=3;
        	case 'U'
        	    x=3;
        	case '-'
        	    x=5;
        	otherwise
        	    x=1;  % NOTE:  Correct to x=5;
		end

	case 4
		switch(a)
	        % full AA alphabet
	        case 'A'
	             x=1;
	        case 'C'    
	            x=2;    
	        case 'G'    
	            x=3;
	        case 'T'
	            x=4;
	        case 'U'
	            x=4;
	        case '-'
	            x=5;
	        otherwise
	            x=1;  % NOTE:  Correct to x=5;
		end
	end
end
end

function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)
% adds pseudocount

    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

    scra = eye(q);

    for i=1:N
        for alpha = 1:q
            for beta = 1:q
               Pij(i,i,alpha,beta) =...
                   (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + ...
                   pseudocount_weight/q*scra(alpha,beta);
            end
        end
    end
end

function C = Compute_C(Pij,Pi,N,q)
% computes correlation matrix

    C=zeros(N*(q-1),N*(q-1));
    for i=1:N
        for j=1:N
            for alpha=1:q-1
                for beta=1:q-1
                     C(mapkey(i,alpha,q),mapkey(j,beta,q)) = ...
                         Pij(i,j,alpha,beta) - Pi(i,alpha)*Pi(j,beta);
                end
            end
        end
    end
end

function A=mapkey(i,alpha,q)
%index mapping to the submatrices
    A = (q-1)*(i-1)+alpha;
end

function W=ReturnW(C,i,j,q)
% extracts coupling matrix for columns i and j

    W = ones(q,q);
    W(1:q-1,1:q-1) = exp( -C(mapkey(i,1:q-1,q),mapkey(j,1:q-1,q)) );

end

function hihj = bp_link(i,j,W,P1,q)
% Adjustment to the correct gauge.

    [mu1, mu2] = compute_mu(i,j,W,P1,q);
    mu1=mu1/mu1(q);
    mu2=mu2/mu2(q);
    hihj=[log(mu1'), log(mu2')];
    return;
end

function [mu1,mu2] = compute_mu(i,j,W,P1,q)
%computes e^hi for each pair using message passing algorithm
    epsilon=1e-4;
    diff =1.0;
    mu1 = ones(1,q)/q;
    mu2 = ones(1,q)/q;
    pi = P1(i,:);
    pj = P1(j,:);

    while ( diff > epsilon )

        scra1 = mu2 * W';
        scra2 = mu1 * W;

        new1 = pi./scra1;
        new1 = new1/sum(new1);

        new2 = pj./scra2;
        new2 = new2/sum(new2);

        diff = max( max( abs( new1-mu1 ), abs( new2-mu2 ) ) );

        mu1 = new1;
        mu2 = new2;

    end
end

function symm=symmetriclocal(Localfield,N,q)
%Copy the information to the bottom part
    symm=Localfield;
    for i=1:N
        for j=(i+1):N
            symm(((j-1)*q+1):j*q,(2*i-1):2*i)=[Localfield(((i-1)*q+1):i*q,2*j),...
                Localfield(((i-1)*q+1):i*q,(2*j-1))];
         end
    end
end

function [hi,sigma]=averagehfield(Pairwisehfield)

%   INPUT:
%       Pairwisehfield  -  q*M x 2*M matrix with the Pairwise h local
%                          fields given by the function DCAparameters.
%   OUTPUT:
%           hi           -  (q-1) x M matrix with the average local fields
%                           by site. WITHOUT GAPS (residues displaced one
%                           place in the numeration i-->(i-1)

    N=size(Pairwisehfield,2)/2;
    q=size(Pairwisehfield,1)/N;
    hi=zeros(q,N);
    sigma=zeros(q,N);
    for i=1:N
        if (i==1)
            hi(1:(q),i)=mean(Pairwisehfield(1:q,3:2:2*N),2);
            sigma(1:(q),i)=std(Pairwisehfield(1:q,3:2:2*N),0,2);
        else if (i==N)
                hi(1:(q),i)=mean(Pairwisehfield(((N-1)*q+1):N*q,1:2:(2*N-3)),2);
                sigma(1:(q),i)=std(Pairwisehfield( ((N-1)*q+1):N*q ,1:2:(2*N-3)),0,2);
            else
                hi(1:(q),i)=mean([Pairwisehfield(( ((i-1)*q+1):((i)*q) ),1:2:(2*(i-1))),...
                    Pairwisehfield( ((i-1)*q+1):i*q ,(2*i+1):2:2*N)],2);
                sigma(1:(q),i)=std([Pairwisehfield(( ((i-1)*q+1):i*q ),1:2:(2*(i-1))),...
                    Pairwisehfield( ((i-1)*q+1):i*q ,(2*i+1):2:2*N)],0,2);
            end
        end
    end
end

function coupling=nicematrix(familycouplinmgs,q)
    n=size(familycouplinmgs,1)/(q-1);
    coupling=zeros(q*n);
    for i=1:n
        for j=1:n
            ii=((i-1)*(q-1)+1):i*(q-1);
            jj=((j-1)*(q-1)+1):j*(q-1);
            newii=((i-1)*q+1):(i*q-1);
            newjj=((j-1)*q+1):(j*q-1);
            coupling(newii,newjj)=familycouplinmgs(ii,jj);
        end
    end
end
