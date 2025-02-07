# print.R
#' print.protein_mods
#' Print method for protein_mods objects
#'
#' @param x A `protein_mods` object to print
#' @param ... other arguments passed to or from 'print'
#'
#' @export
#' @method print protein_mods

print.protein_mods <- function(x, ...)
{
  for(i in 1:dim(x)[1])
  {
    print(x$sequence[[i]])

    if(length(x$mods[[i]]) > 0)
    {
      # identify where each modification is in the sequence
      annotate_at <- ifelse(1:length(x$mods[[i]]) == 1, 4, -1) +
        c(x$mods_at[[i]][1], diff(x$mods_at[[i]]))

      # print out `|` at each modification
      annotate_at |>
        map_chr(~ rep(' ', .x) |>
                  paste(collapse = '') |>
                  paste0('|')) |>
        paste(collapse = '') |>
        paste0('\n') |>
        cat()

      for(j in length(x$mods[[i]]):1)
      {
        # start with `|` at each modification up to the modification we are annotating
        annotate_at[1:j] |>
          map_chr(~ rep(' ', .x) |>
                    paste(collapse = '') |>
                    paste0('|')) |>
          paste(collapse = '') |>
          paste0('\n') |>

          # replace the last `|` with the modification annontation
          str_replace('\\|\n', paste0(x$mods[[i]][j], '\n')) |>

          cat()
      }
    }else{
      cat("  No modifications\n")
    }

    cat("\n")
  }
}

