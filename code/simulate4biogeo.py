#!/usr/bin/env python3

import msprime
import numpy as np
import math
from sys import argv
import subprocess
import gzip

#ts = msprime.simulate(sample_size = 4000, Ne = 0.25, recombination_rate = 10000, mutation_rate = 1000)
#vcffile = open("output.vcf","w")
#ts.write_vcf(vcffile,ploidy = 2)






def write_vcf(tree_sequence, output, ploidy, contig_id):
    """
    Writes a VCF using the sample algorithm as the low level code.
    """
    if tree_sequence.get_sample_size() % ploidy != 0:
        raise ValueError("Sample size must a multiple of ploidy")
    n = tree_sequence.get_sample_size() // ploidy
    sample_names = ["msp_{}".format(j) for j in range(n)]
    last_pos = 0
    positions = []
    for variant in tree_sequence.variants():
        pos = int(round(variant.position))
        if pos <= last_pos:
            pos = last_pos + 1
        positions.append(pos)
        last_pos = pos
    contig_length = int(math.ceil(tree_sequence.get_sequence_length()))
    if len(positions) > 0:
        contig_length = max(positions[-1], contig_length)
    print("##fileformat=VCFv4.2", file=output)
    #print("##source=tskit {}.{}.{}".format(*_tskit.get_tskit_version()), file=output)
    print('##FILTER=<ID=PASS,Description="All filters passed">', file=output)
    print("##contig=<ID={},length={}>".format(contig_id, contig_length), file=output)
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=output)
    print(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", sep="\t", end="", file=output)
    for sample_name in sample_names:
        print("\t", sample_name, sep="", end="", file=output)
    print(file=output)
    for variant in tree_sequence.variants():
        pos = positions[variant.index]
        print(
            contig_id, pos, str(contig_id)+"_"+str(pos), "A", "T", ".", "PASS", ".", "GT",
            sep="\t", end="", file=output)
        for j in range(n):
            genotype = "|".join(
                str(g) for g in
                variant.genotypes[j * ploidy: j * ploidy + ploidy])
            print("\t", genotype, end="", sep="", file=output)
        print(file=output)







def write_vcfgz(tree_sequence, output, ploidy, contig_id):
    """
    Writes a VCF using the sample algorithm as the low level code.
    """
    if tree_sequence.get_sample_size() % ploidy != 0:
        raise ValueError("Sample size must a multiple of ploidy")
    n = tree_sequence.get_sample_size() // ploidy
    sample_names = ["msp_{}".format(j) for j in range(n)]
    last_pos = 0
    positions = []
    for variant in tree_sequence.variants():
        pos = int(round(variant.position))
        if pos <= last_pos:
            pos = last_pos + 1
        positions.append(pos)
        last_pos = pos
    contig_length = int(math.ceil(tree_sequence.get_sequence_length()))
    if len(positions) > 0:
        contig_length = max(positions[-1], contig_length)
    output.write(b"##fileformat=VCFv4.2\n")
    #print("##source=tskit {}.{}.{}".format(*_tskit.get_tskit_version()), file=output)
    output.write(b'##FILTER=<ID=PASS,Description="All filters passed">\n')
    mystr = "##contig=<ID={},length={}>".format(contig_id, contig_length)
    bmystr = bytes(mystr+"\n", 'utf-8')
    output.write(bmystr)
    output.write(b'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    mystr = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    bmystr = bytes(mystr+"", 'utf-8')
    output.write(bmystr)
    for sample_name in sample_names:
        mystr = "\t"+sample_name
        bmystr = bytes(mystr, 'utf-8')
        output.write(bmystr)
    output.write(b'\n')
    #print(file=output)
    for variant in tree_sequence.variants():
        pos = positions[variant.index]
        mystr = str(contig_id)+"\t"+str(pos)+"\t"+str(contig_id)+"_"+str(pos)+"\t"+"A"+"\t"+"T"+"\t"+"."+"\t"+"PASS"+"\t"+"."+"\t"+"GT"
        bmystr = bytes(mystr, 'utf-8')
        output.write(bmystr)
        for j in range(n):
            genotype = "|".join(
                str(g) for g in
                variant.genotypes[j * ploidy: j * ploidy + ploidy])
            bmystr = bytes("\t"+genotype, 'utf-8')
            output.write(bmystr)
        #print(file=output)
        output.write(b'\n')







def symmetric_island_simulation(num_replicates=22,random_seed=None,migration= 1):
    m = migration
    # Allocate the initial sample.
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=800),
        msprime.PopulationConfiguration(sample_size=800),
        msprime.PopulationConfiguration(sample_size=800),
        msprime.PopulationConfiguration(sample_size=800),
        msprime.PopulationConfiguration(sample_size=800)]
    # Now we set up the migration matrix. Since this is a symmetric
    # island model, we have the same rate of migration between all
    # pairs of subpopulations. Diagonal elements must be zero.
    migration_matrix = [
        [0, m, m, m, m],
        [m, 0, m, m, m],
        [m, m, 0, m, m],
        [m, m, m, 0, m],
        [m, m, m, m, 0]]
    # We pass these values to the simulate function, and ask it
    # to run the required number of replicates.
    replicates = msprime.simulate(Ne=0.25, recombination_rate = 10000, mutation_rate = 1000,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        num_replicates=num_replicates,
        random_seed=random_seed)
    # And then iterate over these replicates
    for i, tree_sequence in enumerate(replicates):
        filename = "sym_chromosome_"+str(i)+"_migration_"+str(m)+"_seed_"+str(random_seed)+".vcf.gz"
        #print(filename)
        #vcffile = open(filename,"w")
        #write_vcf(tree_sequence,vcffile,ploidy = 2, contig_id = i)
        #vcffile.close()
        #subprocess.call(["gzip",filename])
        gzipfile = gzip.open(filename, 'wb')
        write_vcfgz(tree_sequence,gzipfile,ploidy = 2,contig_id = i)
        gzipfile.close()



        


#symmetric_island_simulation(22)






def out_of_africa_simulation(num_replicates=22,random_seed=None):
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=1600, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=1600, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=1600, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    replicates = msprime.simulate(recombination_rate = 1.25, mutation_rate = 0.125,  #  0.01 --> ca 8000 SNPs
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        num_replicates=num_replicates,
        random_seed=random_seed)
    # And then iterate over these replicates
    for i, tree_sequence in enumerate(replicates):
        filename = "ooa_chromosome_"+str(i+1)+"_seed_"+str(random_seed)+".vcf.gz"
        #vcffile = open(filename,"w")
        #write_vcf(tree_sequence,vcffile,ploidy = 2,contig_id = i)
        #vcffile.close()
        #subprocess.call(["gzip",filename])
        gzipfile = gzip.open(filename, 'wb')
        write_vcfgz(tree_sequence,gzipfile,ploidy = 2,contig_id = i)
        gzipfile.close()
        
#out_of_africa_simulation(1)








if argv[1]=="si":
    migration = float(argv[2])
    seed = int(argv[3])
    num_chromosomes = int(argv[4])
    #print("simulating symetric islands\n")
    symmetric_island_simulation(num_chromosomes,seed,migration)


if argv[1]=="ooa":
    seed = int(argv[2])
    num_chromosomes = int(argv[3])
    print("simulating ooa")
    out_of_africa_simulation(num_chromosomes,seed)





#how to call symmetric_island_simulation from within R
#globalcommand = paste("python simulate4biogeo si ", migration , " ", seed ," ", chromosomes ," ", sep="")
#globalcommand = "python simulate4biogeo si "
#commands = c()
#seeds = sample(5000,length(migrationrates))
#for(i in 1:length(migrationrates)){
    #commands = c(commands,paste(globalcommand, migrationrates[i],seeds[i],"3", sep= " "))
#}

#f<-function(commandline){
    #system(commandline, ignore.stdout=TRUE,ignore.stderr=TRUE,wait=FALSE)
#}
#library(parallel)
#mclapply(commands,f,mc.cores=4)




#how to call out_of_africa_simulation from within R
#globalcommand = "python simulate4biogeo ooa "
#commands = c()
#for(seed in 1:3){
    #commands = c(commands,paste(globalcommand, seed, sep= " "))
#}

#f<-function(commandline){
    #system(commandline, ignore.stdout=TRUE,ignore.stderr=TRUE,wait=FALSE)
#}
#library(parallel)
#mclapply(commands,f,mc.cores=4)
