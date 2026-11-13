# Capture R coloc.abf reference numbers for tests/fixtures/coloc_abf_r_reference.json.
#
# The R `coloc` package is GPL and ClawBio is MIT, so it is never installed,
# imported or vendored into this repository, and this script is NOT run by the
# test suite. It is here so that anyone can re-derive the reference numbers the
# parity test asserts against, rather than taking them on trust.
#
# Run it once, outside this repository's environment, paste the printed vector
# into the fixture under "r_reference", and set "present": true.
#
#   Rscript capture_r_reference.R         # needs install.packages("coloc")
#
# The inputs below must stay identical to the fixture's "inputs" block (the same
# snp, beta, se, sdY and priors), or the comparison means nothing.

library(coloc)

snp  <- c("1_109000000_A_G","1_109000500_A_G","1_109001000_A_G","1_109001500_A_G",
          "1_109002000_A_G","1_109002500_A_G","1_109003000_A_G","1_109003500_A_G",
          "1_109004000_A_G","1_109004500_A_G","1_109005000_A_G","1_109005500_A_G")
beta1 <- c(0.01,-0.02,0.18,0.03,-0.01,0.22,0.05,-0.03,0.02,0.01,-0.02,0.04)
se1   <- c(0.04,0.04,0.035,0.04,0.05,0.03,0.045,0.04,0.04,0.05,0.04,0.04)
beta2 <- c(0.00,-0.01,0.15,0.02,0.01,0.20,0.03,-0.02,0.01,0.02,-0.01,0.03)
se2   <- c(0.05,0.05,0.04,0.05,0.05,0.035,0.05,0.05,0.05,0.05,0.05,0.05)

D1 <- list(beta = beta1, varbeta = se1^2, snp = snp, type = "quant", sdY = 1.0)
D2 <- list(beta = beta2, varbeta = se2^2, snp = snp, type = "quant", sdY = 1.0)

res <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

cat("Paste these into r_reference, and set present = true:\n")
print(res$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")])
