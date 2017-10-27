context("test3sm")

test_that("test3sm() output", {
  dipper = system.file("extdata", "ed.inp", package = "R2ucare")
  dipper = read_inp(dipper,group.df=data.frame(sex=c('Male','Female')))
  dip.hist = dipper$encounter_histories
  dip.freq = dipper$sample_size
  dip.group = dipper$groups
  mask = (dip.group == 'Female')
  dip.fem.hist = dip.hist[mask,]
  dip.fem.freq = dip.freq[mask]
  res.females = test3sm(dip.fem.hist, dip.fem.freq)
  expect_output(str(res.females), "List of 2")
  expect_equivalent(res.females$test3sm, c(2.041,3,0.564))
})

