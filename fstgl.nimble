# Package metadata
version = "0.1.0"
author = "Andy Jasonowicz"
description = "A program to estimate fst from genotype likelihoods"

# Dependencies
requires "nim >= 1.0.0"
requires "hts >= 0.3.22"
requires "zip"
requires "therapist"


bin = @["src/glUtils"]
