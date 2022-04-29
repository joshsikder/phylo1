import msprime, newick
from IPython.display import SVG, display

def simpleTree(numSamp): 
    tree = msprime.sim_ancestry(numSamp)
    return tree

def simplePopulations(N, samples):
    demography = msprime.Demography.stepping_stone_model([100]*N, migration_rate=0.1, boundaries=True)
    tree = msprime.sim_ancestry({0: samples, N-1: samples}, demography=demography)
    return tree

def speciesTree(ploidy, seqlength):
    demography = msprime.Demography.from_species_tree("(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)",time_units="myr",initial_size=10**4,generation_time=20)
    tree = msprime.sim_ancestry({"gorilla": 2, "human": 4}, demography=demography, ploidy=ploidy, sequence_length=seqlength)
    return tree

def firstSim(seqlength):
    time_units = 1000/25 # kya -> generations
    demo1 = msprime.Demography()
    N=10**4
    demo1.add_population(name="Dsantomea", initial_size=N)
    demo1.add_population(name="Dyakuba", initial_size=N)
    demo1.add_population(name="Derecta", initial_size=N)
    demo1.add_population_split(time=5, derived=["Dsantomea","Dyakuba"],ancestral="Derecta")
    demo1.add_mass_migration(time=50*time_units, source="Dsantomea", dest="Derecta", proportion=.1)
    treeseq = msprime.sim_ancestry(
        recombination_rate=1e-8,
        sequence_length=seqlength,
        samples=[
            msprime.SampleSet(1, ploidy=1, population="Dsantomea"),
            msprime.SampleSet(1, ploidy=1, population="Dyakuba"),
            msprime.SampleSet(1, ploidy=1, population="Derecta")
        ],
        demography=demo1,
        record_migrations=True,
        random_seed=123456
    )
    print(f"Simulation of {sequence_length/10**3}Kb run, using record_migrations=True")
    print(
        "NB: time diff from Neanderthal split to admixture event is",
        f"{300 * time_units - 50 * time_units:.0f} gens",
        f"({(300 * time_units - 50 * time_units) / 2 / Ne} coalescence units)"
    )
    return ts

#ts = speciesTree(1, 10000)
ts = firstSim(30 * 10**3)
ts.write_fasta('/nas/longleaf/home/jsikder/schriderlab/test.fasta')
with open('/nas/longleaf/home/jsikder/schriderlab/tree.draw') as outfile:
    outfile.write(ts.draw_text())
with open('/nas/longleaf/home/jsikder/schriderlab/tree.stats') as outfile:
    outfile.write(ts)