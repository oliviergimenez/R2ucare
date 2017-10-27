context("test2cl")

test_that("test2cl() output", {
  dipper = system.file("extdata", "ed.inp", package = "R2ucare")
  dipper = read_inp(dipper,group.df=data.frame(sex=c('Male','Female')))
  dip.hist = dipper$encounter_histories
  dip.freq = dipper$sample_size
  dip.group = dipper$groups
  mask = (dip.group == 'Male')
  dip.mal.hist = dip.hist[mask,]
  dip.mal.freq = dip.freq[mask]
  X = dip.mal.hist
  freq = dip.mal.freq
  res.males = test2cl(X,freq)
  expect_output(str(res.males), "List of 2")
  expect_equivalent(res.males$test2cl, c(0,0,1))
})

