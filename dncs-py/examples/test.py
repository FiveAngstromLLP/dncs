import dncs

# Create a polymer object
polymer = dncs.Polymer("YGGFM")
# Energy of the polymer
print("Energy = ", polymer.getEnergy())

# Create a sampler object
sample = dncs.SobolSampler(polymer, no_of_samples=100, sidechain=False)


# Sample the polymer
sample.toPDB("sample.pdb")
# sample.toPDBFiles("sample")
