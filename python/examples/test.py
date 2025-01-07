import dncs

# Create a polymer object
polymer = dncs.Polymer("YGGFM", "amberfb15.xml")
# Energy of the polymer
# print("Energy = ", polymer.getEnergy())

# Create a sampler object
sample = dncs.SobolSampler(polymer, 100, 4, "Result/Sample")
