import dncs

# Create a polymer object
polymer = dncs.Polymer("YGGFM")
# Energy of the polymer
print("Energy = ", polymer.getEnergy())

# Create a sampler object
sample = dncs.SobolSampler(polymer, no_of_samples=15, sidechain=False)

# Create a Minimizer object
# minimizer = dncs.Minimizer(sample)
# # Minimize the polymer
# minimizer.minimize()


# Sample the polymer
sample.toPDB("sample.pdb")
sample.toPDBFiles("sample")
