context("test2ct")

test_that("test2ct() output", {
  dipper = system.file("extdata", "ed.inp", package = "R2ucare")
  dipper = read_inp(dipper,group.df=data.frame(sex=c('Male','Female')))
  dip.hist = dipper$encounter_histories
  dip.freq = dipper$sample_size
  dip.group = dipper$groups
  mask = (dip.group == 'Female')
  dip.fem.hist = dip.hist[mask,]
  dip.fem.freq = dip.freq[mask]
  X = dip.fem.hist
  freq = dip.fem.freq
  res.females = test2ct(X,freq)
  expect_output(str(res.females), "List of 2")
  expect_equivalent(res.females$test2ct[1:3], c(3.250,4,0.517))
})

