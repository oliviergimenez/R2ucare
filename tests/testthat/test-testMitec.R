context("testMitec")

test_that("testMitec() output", {
  geese = system.file("extdata", "geese.inp", package = "R2ucare")
  geese = read_inp(geese)
  geese.hist = geese$encounter_histories
  geese.freq = geese$sample_size
  skip_on_cran()
  res = testMitec(geese.hist,geese.freq)
  expect_output(str(res), "List of 2")
  expect_equivalent(res$testMitec[2:3], c(27,0))
})

