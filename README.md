# DCAscapes
Codebase for Direct Coupling Analysis applied to nucleotide recognition

This package contains codes specific for DCA-scape analysis:

	
	1. masterDCAparameters_length - estimating the mean couplings and local field parameters from the input nucleotide fasta data for each
 	of the four nucleotides gauged states. 
	2.Fastahamiltonian - Calculating the Hamiltonian score for any 20mer sequences list base on the average couplings and local fields. 
	Lambda genome sequence fasta file (lamdba_genome_seq.fasta) is provided as an example for the Hamiltonian score calculation. 

All of these codes are based on MATLAB. Below is one example for the basic analysis of lambda SEQRS data. 

[c_average,h_average] = masterDCAparameters_length('RPI_033_lambda.fasta',20);
[c_average_GST,h_average_GST] = masterDCAparameters_length('RPI_033_GST_lambda.fasta',20);

Hamiltonian_score_033=Fastahamiltonian('lambda_genome_seq.fasta',c_average_033,h_average_033,2,0,2);
Hamiltonian_score_gst=Fastahamiltonian('lambda_genome_seq.fasta',c_average_031_GST,h_average_031_GST,2,0,2);

Hamiltonian_score_033_gst=Hamiltonian_score_033-Hamiltonian_score_gst;

Any publication resulting from applications of DCA-scape should cite:

Q Zhou, N Kunder, José Alberto De la Paz, AE. Lasley, VD.Bhat, 
F Morcos, ZT. Campbell (2018),Global pairwise RNA interaction 
landscapes reveal corefeatures of protein recognition.
