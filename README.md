# DCAscapes
Codebase for Direct Coupling Analysis applied to nucleotide recognition

This package contains codes specific for DCA-scape analysis:

	
	1.DCAparameters    - estimating the couplings and local fields parameters from the input SEQRS fasta data. There are four subtypes of DCAparameters 
	estimation codes based on the nucleotide used as a reference to reduce the freedom form independent constraints.
	2.average_couplings_localfields    - Calculating the mean of the couplings and local fields for each of the four nucleotides gauged states. 
	3.Fastahamiltonian    - Calculating the Hamiltonian score for any 20mer sequences list base on the average couplings and local fields. 
	Lambda genome sequence fasta file (lamdba_genome_seq.fasta) is provided as an example for the Hamiltonian score calculation. 

All of these codes are based on MATLAB. Below is one example for the basic analysis of lambda SEQRS data. 

[h4_033,familycouplings4_033,align4_033]=DCAparameters_4('RPI_033_lambda.fasta',2);
[h3_033,familycouplings3_033,align3_033]=DCAparameters_3('RPI_033_lambda.fasta',2);
[h2_033,familycouplings2_033,align2_033]=DCAparameters_2('RPI_033_lambda.fasta',2);
[h1_033,familycouplings1_033,align1_033]=DCAparameters_1('RPI_033_lambda.fasta',2);


[h4_031_GST,familycouplings4_031_GST,align4_031_GST]=DCAparameters_4('RPI_031_GST_lambda_p22.fasta',2);
[h3_031_GST,familycouplings3_031_GST,align3_031_GST]=DCAparameters_3('RPI_031_GST_lambda_p22.fasta',2);
[h2_031_GST,familycouplings2_031_GST,align2_031_GST]=DCAparameters_2('RPI_031_GST_lambda_p22.fasta',2);
[h1_031_GST,familycouplings1_031_GST,align1_031_GST]=DCAparameters_1('RPI_031_GST_lambda_p22.fasta',2);

[c_average_033,h_average_033] = average_couplings_localfields(familycouplings1_033,familycouplings2_033,familycouplings3_033,familycouplings4_033,h1_033,h2_033,h3_033,h4_033);
[c_average_031_GST,h_average_031_GST] = average_couplings_localfileds(familycouplings1_031_GST,familycouplings2_031_GST,familycouplings3_031_GST,familycouplings4_031_GST,h1_031_GST,h2_031_GST,h3_031_GST,h4_031_GST);

Hamiltonian_score_033=Fastahamiltonian('lambda_genome_seq.fasta',c_average_033,h_average_033,2,0,2);
Hamiltonian_score_gst=Fastahamiltonian('lambda_genome_seq.fasta',c_average_031_GST,h_average_031_GST,2,0,2);

Hamiltonian_score_033_gst=Hamiltonian_score_033-Hamiltonian_score_gst;

Any publication resulting from applications of DCA-scape should cite:

Q Zhou, N Kunder, José Alberto De la Paz, AE. Lasley, VD.Bhat, 
F Morcos, ZT. Campbell (2018),Global pairwise RNA interaction 
landscapes reveal corefeatures of protein recognition.
