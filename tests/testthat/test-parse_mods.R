test_that("multiplication works", {
  data(MSfrags)

  # check with one string
  test1 <- parse_mods(MSfrags[1])

  expect_s3_class(test1, "protein_mods")
})
