rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_range100.html",
                  output_dir = "doc",
                  params = list(range0 = 100),
                  clean = TRUE,
                  envir = new.env()
)

rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_range500.html",
                  output_dir = "doc",
                  params = list(range0 = 500),
                  clean = TRUE,
                  envir = new.env()
)

rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_range1000.html",
                  output_dir = "doc",
                  params = list(range0 = 1000),
                  clean = TRUE,
                  envir = new.env()
)

rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_range5000.html",
                  output_dir = "doc",
                  params = list(range0 = 5000),
                  clean = TRUE,
                  envir = new.env()
)