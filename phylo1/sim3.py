import msprime

def ancestrySim(seqlength, randomseed):
    # time_units = 4 # kya -> generations
    speciesTree = msprime.Demography()
    N=1000000 #Ne
    time_units = 4*N
    speciesTree.add_population(name="Dsantomea", initial_size=N)
    speciesTree.add_population(name="Dyakuba", initial_size=N)
    speciesTree.add_population(name="Derecta", initial_size=N)
    speciesTree.add_population(name="santomea_yakuba", initial_size=N)
    speciesTree.add_population(name="ancestor", initial_size=N)
    speciesTree.add_mass_migration(time=0.1*time_units, source="Dyakuba", dest="Derecta", proportion=.1)
    speciesTree.add_population_split(time=0.5*time_units, derived=["Dsantomea","Dyakuba"],ancestral="santomea_yakuba")
    speciesTree.add_population_split(time=1*time_units, derived=["santomea_yakuba","Derecta"],ancestral="ancestor")
    
    treeseq = msprime.sim_ancestry(
        recombination_rate=1e-8,
        # recombination_rate = 0,
        sequence_length=seqlength,
        samples=[
            msprime.SampleSet(1, ploidy=1, population="Dsantomea"),
            msprime.SampleSet(1, ploidy=1, population="Dyakuba"),
            msprime.SampleSet(1, ploidy=1, population="Derecta")
        ],
        demography=speciesTree,
        record_migrations=True,
        random_seed=randomseed
    )
    return treeseq

ts = ancestrySim(100, 123456)
model = msprime.MatrixMutationModel(
    ["A","T","G","C"],
    root_distribution = [0.25, 0.25, 0.25, 0.25],
    transition_matrix = [[0.0, 1/3, 1/3, 1/3],
                         [1/3, 0.0, 1/3, 1/3],
                         [1/3, 1/3, 0.0, 1/3],
                         [1/3, 1/3, 1/3, 0.0]]
)
msprime.sim_mutations(ts, rate=2, model=model)
ts.write_fasta('/nas/longleaf/home/jsikder/schriderlab/ignored/phylo1/bintro.fasta')
def write_output():
    with open('/nas/longleaf/home/jsikder/schriderlab/ignored/phylo1/bintro.draw', "w+") as outfile:
        outfile.write(ts.draw_text())
write_output()
#print(ts.draw_text())