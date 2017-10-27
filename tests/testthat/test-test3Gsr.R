context("test3Gsr")

test_that("test3Gsr() output", {
  geese = system.file("extdata", "geese.inp", package = "R2ucare")
  geese = read_inp(geese)
  geese.hist = geese$encounter_histories
  geese.freq = geese$sample_size
  res = test3Gsr(geese.hist,geese.freq)
  expect_output(str(res), "List of 2")
  expect_equivalent(res$test3Gsr, c(117.753,12,0))
})

