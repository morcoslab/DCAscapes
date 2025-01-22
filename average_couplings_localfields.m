function [c_average,h_average] = average_couplings_localfields(familycouplings1,familycouplings2,familycouplings3,familycouplings4,h1,h2,h3,h4)
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



% transfer clouplings 1 to 4
c_4_1=zeros(80);
for (i=[1:4:80])
    for (j=[1:4:80])
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
c_4_2=zeros(80);
for (i=[1:4:80])
    for (j=[1:4:80])
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
c_4_3=zeros(80);
for (i=[1:4:80])
    for (j=[1:4:80])
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



h_4_1=zeros(4,20);
%%% transfer 1 to 4
for i=1:20
    h_4_1(1,i)=h1(4,i);
    h_4_1(2,i)=h1(2,i);
    h_4_1(3,i)=h1(3,i);
    h_4_1(4,i)=h1(1,i);
end

h_4_2=zeros(4,20);
%%% transfer 2 to 4
for i=1:20
    h_4_2(1,i)=h2(1,i);
    h_4_2(2,i)=h2(4,i);
    h_4_2(3,i)=h2(3,i);
    h_4_2(4,i)=h2(1,i);
end

h_4_3=zeros(4,20);
%%% transfer 3 to 4
for i=1:20
    h_4_3(1,i)=h3(1,i);
    h_4_3(2,i)=h3(2,i);
    h_4_3(3,i)=h3(4,i);
    h_4_3(4,i)=h3(3,i);
end

h_average=(h_4_1+h_4_2+h_4_3+h4)/4;
