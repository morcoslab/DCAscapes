function H=Fastahamiltonian(inputfile,couplings,localfields,Htype,N1,stype)
% Direct Coupling Analysis (DCA) Scapes - Fastahamiltonian scores calculation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	General Hamiltonian score H
%	INPUT:
%		inputfile   - file containing the FASTA alignment.
%		couplings   - Corresponding coupling matrix as returned
%			      from the DCAparameters function of the DCA
%			      framework.
%		localfields - Corresponding h fields matrix as returned
%			      from the DCAparameters function of the DCA
%			      framework.
%		Htype 	    - Hamiltonian type: type 1 sums couplings over
%			      pairs across two species while type 2
%			      evaluates the complete Hamiltonian for the
%			      complete sequence as in the traditional Potts
%			      model.
%		N1 	    - Lenght of the first species. This number is
%			      ignored in the Htype 2 case.
%		stype 	    - species type: 1 for proteins
%			   		    2 for RNA and DNA
%
%	OUTPUT:
%		H 	    - Vector of the Hamiltonians of each sequences.
%
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
%  Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    [N,M,q,Sequences] = return_alignment(inputfile, stype);
    H=zeros(M,1);

    %Local Fields
    for seq=1:M
        for res=1:N
            H(seq)=localfields((Sequences(seq,res)),res)+H(seq);
        end
    end

    %Couplings
    switch Htype
        case 1
            %First residue in the first domain
            %Compute across 2 different species
              for seq=1:M
                  for res=1:N1
                      for pair=(N1+1):N
                          iindex=q*(res-1)+Sequences(seq,res);
                          jindex=q*(pair-1)+Sequences(seq,pair);
                          H(seq)=couplings(iindex,jindex)+H(seq);
                      end
                  end
              end
        case 2
            %Complete Hamiltonian
            for seq=1:M
                  for res=1:N
                      if (res<N)
                          for pair=(res+1):N
                              iindex=q*(res-1)+Sequences(seq,res);
                              jindex=q*(pair-1)+Sequences(seq,pair);
                              H(seq)=couplings(iindex,jindex)+H(seq);
                          end
                      end
                  end
            end
    end
    H=-H;
end

function [N,M,q,Z] = return_alignment(inputfile,type)
% reads alignment from inputfile, removes inserts and converts into numbers

    align_full = fastaread(inputfile);
    M = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Z = zeros(M,N);

    for i=1:M
        counter = 0;
        for j=1:length(ind)
            if( ind(j) )
                counter = counter + 1;
                Z(i,counter)=letter2number( align_full(i).Sequence(j),type );
            end
        end
    end
    q = max(max(Z));
end

function x=letter2number(a,type)
	switch (type)
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
            x=1;
    end
    end
end
