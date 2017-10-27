context("test3sr")

test_that("test3sr() output", {
  dipper = system.file("extdata", "ed.inp", package = "R2ucare")
  dipper = read_inp(dipper,group.df=data.frame(sex=c('Male','Female')))
  dip.hist = dipper$encounter_histories
  dip.freq = dipper$sample_size
  dip.group = dipper$groups
  mask = (dip.group == 'Male')
  dip.mal.hist = dip.hist[mask,]
  dip.mal.freq = dip.freq[mask]
  res.males = test3sr(dip.mal.hist, dip.mal.freq)
  expect_output(str(res.males), "List of 2")
  expect_equivalent(res.males$test3sr[1:2], c(6.778,5))
})

