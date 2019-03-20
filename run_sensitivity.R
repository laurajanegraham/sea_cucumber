rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_uninf",
                  output_dir = "doc",
                  clean = TRUE,
                  envir = new.env()
                  )

rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_uninf",
                  output_dir = "doc",
                  params = list(pcprior = TRUE),
                  clean = TRUE,
                  envir = new.env()
)

rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_uninf",
                  output_dir = "doc",
                  params = list(pcprior = TRUE, range0 = 50),
                  clean = TRUE,
                  envir = new.env()
)

rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_uninf",
                  output_dir = "doc",
                  params = list(pcprior = TRUE, range0 = 500),
                  clean = TRUE,
                  envir = new.env()
)

rmarkdown::render(input = "sea_cucumber_spatial_hurdle.Rmd", 
                  output_format = "html_document",
                  output_file = "sea_cucumber_uninf",
                  output_dir = "doc",
                  params = list(pcprior = TRUE, range0 = 5000),
                  clean = TRUE,
                  envir = new.env()
)