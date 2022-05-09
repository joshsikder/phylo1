import msprime, newick
from IPython.display import SVG, display

def firstSim(seqlength):
    time_units = 1000/25 # kya -> generations
    demo1 = msprime.Demography()
    N=1000
    demo1.add_population(name="A", initial_size=N)
    demo1.add_population(name="B", initial_size=N)
    demo1.add_population(name="C", initial_size=N)
    demo1.add_population(name="D", initial_size=N)
    demo1.add_population(name="E", initial_size=N)
    demo1.add_population(name="F", initial_size=N)
    demo1.add_population_split(time=50, derived=["ABC","DEF"],ancestral="ABCDEF")
    demo1.add_population_split(time=30, derived=["AB","C"],ancestral="ABC")
    demo1.add_population_split(time=25, derived=["D","EF"],ancestral="DEF")
    demo1.add_population_split(time=10, derived=["A","B"],ancestral="AB")
    demo1.add_population_split(time=10, derived=["E","F"],ancestral="EF")
    demo1.add_mass_migration(time=5, source="B", dest="C", proportion=.02)
    demo1.add_mass_migration(time=25, source="ABC", dest="D", proportion=.07)
    demo1.add_mass_migration(time=15, source="D", dest="C", proportion=.1)

    treeseq = msprime.sim_ancestry(
        recombination_rate=1e-16,
        sequence_length=seqlength,
        samples=[
            msprime.SampleSet(1, ploidy=1, population="A"),
            msprime.SampleSet(1, ploidy=1, population="B"),
            msprime.SampleSet(1, ploidy=1, population="C"),
            msprime.SampleSet(1, ploidy=1, population="D"),
            msprime.SampleSet(1, ploidy=1, population="E"),
            msprime.SampleSet(1, ploidy=1, population="F")
        ],
        demography=demo1,
        record_migrations=True,
        random_seed=123456
    )
    print(f"Simulation of {seqlength/10**3}Kb run, using record_migrations=True")
    print(
        "NB: time diff from Neanderthal split to admixture event is",
        f"{300 * time_units - 50 * time_units:.0f} gens",
        f"({(300 * time_units - 50 * time_units) / 2 / N} coalescence units)"
    )
    return treeseq

ts = firstSim(1 * 10**3)
ts.write_fasta('/nas/longleaf/home/jsikder/schriderlab/ignored/phylo1/multiintro.fasta')
with open('/nas/longleaf/home/jsikder/schriderlab/ignored/phylo1/multiintro.draw') as outfile:
    outfile.write(ts.draw_text())
with open('/nas/longleaf/home/jsikder/schriderlab/ignored/phylo1/multiintro.stats') as outfile:
    outfile.write(ts)