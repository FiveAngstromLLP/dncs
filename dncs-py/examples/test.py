import dncs

polymer = dncs.Polymer("YGGFM")

sample = dncs.SobolSampler(polymer, 10)


sample.toPDB("sample.pdb")
