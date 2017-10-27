context("test3Gwbwa")

test_that("test3Gwbwa() output", {
  geese = system.file("extdata", "geese.inp", package = "R2ucare")
  geese = read_inp(geese)
  geese.hist = geese$encounter_histories
  geese.freq = geese$sample_size
  res = test3Gwbwa(geese.hist,geese.freq)
  expect_output(str(res), "List of 2")
  expect_equivalent(res$test3Gwbwa, c(472.855,20,0))
})

