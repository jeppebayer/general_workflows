//Number of population samples (demes)
2
//Population effective sizes (number of genes)
NPOP1
NPOP2
//Samples sizes and samples age 
100
100
//Growth rates	: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 mig_rec
mig_rec 0
//Migration matrix 1
0 mig_int
mig_int 0
//Migration matrix 2
0 mig_anc
mig_anc 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
3  historical event
600 0 0 0 1 0 1
XXX 0 0 0 1 0 2
TDIV 1 0 1 RESIZE0 0 2
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   5e-9

