#' @export
print.cqc_cf_list <- function(x, ...) {
  cat("cytoqc data: \n")
  cat(length(x), "samples \n")
}

#' @export
print.cqc_gs <- function(x, ...) {
  cat("cytoqc data: \n")
  cat("List of GatingHierarchy objects with", length(x), "samples \n")
}

#' @importFrom knitr knit_print
#' @export
knit_print.cqc_reference <- function(x, ...) {
  x %>%
    kable() %>%
    kable_styling("bordered", full_width = F, position = "left", font_size = 12) %>%
    row_spec(0, background = "gray", color = "black") %>%
    knit_print()
}

#' @export
print.cqc_match_result <- function(x, ...) {
  print(as_tibble(x))
}

#' @export
print.cqc_match_result_and_solution <- function(x, ...) {
  print(format(x, ...))
}


#' @export
knit_print.cqc_match_result_and_solution <- function(x, ...) {
  match_result_to_dt(x,...) %>% knit_print

}

#' color encoding table for the match report
#' @noRd
match_result_color_tbl <- function(x, ...) {

  df_color <- format(x, ...)#check mark preserved version of coloring
  ncol <- ncol(df_color)
  df_color[] <- lapply(df_color, as.character)
  #start color encoding
  ref_vec <- df_color[,"Ref"]
  nref <- sum(ref_vec!="")
  df_color[1:nref,"Ref"] <- "black"
  for(i in seq_len(ncol)[-(1:2)])
  {
    vec <- df_color[1:nref, i]
    vec[is.na(vec)] <- ""
    #color the approximately matched
    vec[vec != "" & vec != "\u2713"] <- "green"
    #color the extact matched
    vec[vec == "\u2713"] <- "gray"
    #color the unmatched
    vec[vec==""] <- "red"
    df_color[1:nref, i] <- vec
  }

  ridx.Unmatched <- which(df_color[,1]=="Unmatched")
  if(length(ridx.Unmatched) > 0)
  {
    vec <- df_color[ridx.Unmatched, ]
    row <- vec[, -c(1:2)]
    row[row != ""] <- "red"
    vec[, -c(1:2)] <- row
    vec[1] <- "black"
    df_color[ridx.Unmatched, ] <- vec
  }

  ridx.rm <- which(df_color[,1]=="To Delete")
  if(length(ridx.rm) > 0)
  {
    vec <- df_color[ridx.rm, ]
    row <- vec[, -c(1:2)]
    row[row != ""] <- "brown"
    vec[, -c(1:2)] <- row
    vec[1] <- "black"
    df_color[ridx.rm, ] <- vec
  }

  df_color[df_color==""] <-  "white"
  df_color
}
#' This can not be called at rstudio console since knit_print requires PhantomJS when run interactively
#' @importFrom DT formatStyle styleEqual datatable
#' @noRd
match_result_to_dt <- function(x, ...) {
  df <- format(x, show_check_mark = FALSE, ...)
  df[is.na(df)] <- "N/A" #display NA in dt
  df_color <- match_result_color_tbl(x, ...)
  ncol <- ncol(df_color)

  #append df_color
  df <- cbind(df, df_color)
  col_ref_col_idx <- (ncol+1):(ncol+ncol)
  #apply the color map
  colors <- unique(unlist(df_color))

  datatable(df
            , class = "compact"
            , colnames = c("", colnames(df)[-1])
            , filter = "none"
            , options = list(
     columnDefs = list(list(targets = col_ref_col_idx, visible = FALSE)) #hide df_color cols
     # , pageLength = 12
     , paging = FALSE
     , searching = FALSE
     , info = FALSE
     , ordering = FALSE
     , dom = 't'
   ))%>%
          formatStyle(1:ncol, col_ref_col_idx
                      , color = styleEqual(colors, colors)
                      # , fontSize = "80%"
          )


}
#' @export
#' @importFrom purrr map_dfc
format.cqc_match_result_and_solution <- function(x, show_check_mark = TRUE, ...) {
  #combine match res and recommened solution to present it as wide format for easy viewing the data
  ref <- x[["ref"]]
  match_result <- x[["match_result"]]
  gids <- names(x[["match_result"]])
  tbl <- map_dfc(gids, function(gid) {
      df <- filter(x[["solution"]], group_id == gid)
      this_res <- match_result[[gid]]
      #init column
      if(show_check_mark)
        col_to_show <- rep("\u2713", length(ref))
      else
        col_to_show <- ref
      if(!is.null(this_res))
      {
        missing <- this_res[["missing"]]
        unknown <- this_res[["unknown"]]
         #drop the recommended deletion
        df1 <- filter(df, !is.na(to))
        #drop the recommended insertion
        df1 <- filter(df1, !is.na(from))
        #fill the refs with the recommended edit
        matched.ref <- df1[["to"]]
        matched.target <- df1[["from"]]
        idx <- match(matched.ref, ref)
        col_to_show[idx] <- matched.target
        # browser()
        #fill the unmatched refs
        unmatched.ref <- missing[!missing%in%matched.ref]
        col_to_show[match(unmatched.ref, ref)] <- NA

        ## the redundant item
        torm <- filter(df, is.na(to))[["from"]]

        #insertion item (may overwrite unmatched entry of NA)
        insert.ref <- filter(df, is.na(from))[["to"]]
        col_to_show[match(insert.ref, ref)] <- "<TO INSERT>"

        # the unmatched target
        if(length(matched.target)>0)
          unmatched.target <- unknown[-match(matched.target, unknown)]#exclude the matched ones
        else
          unmatched.target <- unknown

        unmatched.target <- unmatched.target[!unmatched.target %in% torm]#exclude the deletion entries

        #append unmatched.target
        if(length(unmatched.target) == 0)
          unmatched.target <- ""
        else
          unmatched.target <- paste(unmatched.target, collapse = ",")
        col_to_show <- c(col_to_show, unmatched.target)
        #append deletions
        if(length(torm) == 0)
          torm <- ""
        else
          torm <- paste(torm, collapse = ",")
        col_to_show <- c(col_to_show, torm)
      }else
        col_to_show <- c(col_to_show, "", "")

      col_to_show <- list(col_to_show)
      names(col_to_show) <- gid
      col_to_show
    })
  tbl <- cbind(c(ref, "", ""), tbl)
  tbl <- cbind(c(rep("", length(ref)), "Unmatched", "To Delete"), tbl)
  colnames(tbl)[1:2] <- c("", "Ref")

   #rm last two rows if they are all empty
  ridx <- nrow(tbl) - 1
  unmatched <- tbl[ridx, -1]
  if(isTRUE(all(unmatched == "")))
    tbl <- tbl[-ridx,]

  ridx <- nrow(tbl)
  torm <- tbl[ridx,-1]
  if(isTRUE(all(torm == "")))
    tbl <- tbl[-ridx,]

  tbl
}
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling row_spec
#' @export
knit_print.cqc_match_result <- function(x, ...) {
  if (length(x) == 0) {
    res <- kable(data.frame(object = "All passed"), col.names = NULL) %>% row_spec(1, color = "green")
  } else {
    res <- kable(as_tibble(x)) %>%
      row_spec(0, background = "#9ebcda", color = "black")
  }

  res %>%
    kable_styling("bordered", full_width = F, position = "left", font_size = 12) %>%
    column_spec(1, bold = TRUE) %>%
    knit_print()
}

#' @importFrom htmltools htmlEscape
#' @importFrom kableExtra collapse_rows cell_spec column_spec
#' @importFrom dplyr %>% select distinct mutate_if
#' @importFrom tidyr unite
#' @export
knit_print.cqc_solution <- function(x, itemize = FALSE, ...) {
  if (!itemize) {
    x <- x %>%
      select(-1) %>%
      distinct()
  }

  x <- x %>% mutate(from = htmlEscape(from), to = htmlEscape(to)) %>%
    mutate("Proposed change" = ifelse(is.na(to), cell_spec(from, strikeout = TRUE), paste(from, to, sep = " --> "))) %>%
    select(-c(from, to)) %>%
    kable(escape = F) %>%
    kable_styling(c("bordered", "condensed"), full_width = F, position = "left", font_size = 12) %>%
    collapse_rows(columns = 1, "top") %>%
    row_spec(0, background = "#e5f5e0", color = "black")
  if (itemize) {
    x <- x %>% column_spec(1, bold = TRUE)
  }
  knit_print(x)
}
#' @importFrom dplyr summarise
collapse_params <- function(x, ...) {
  class_names <- class(x)
  type <- sub("cqc_check_", "", class_names[2])
  if (type != "panel") {
    type <- as.symbol(type)
    x <- group_by(x, group_id, nObject) %>%
      summarise(!!type := paste(!!type, collapse = ", ")) %>%
      arrange(desc(nObject))
    class(x) <- class_names
  }


  x
}

#' @export
print.cqc_check <- function(x, ...){
   print(summary(x), ...)
}

#' @noRd
#' @importFrom knitr is_latex_output is_html_output
#' @export
knit_print.cqc_cluster_panel <- function(x, ...){
  # browser()
  check_res <- attr(x, "check_res")
  value <- "marker"

  summarized <- check_res %>% summary

  df <- summarized %>%
    select(-nObject) %>%
    arrange(group_id) %>%
    spread(group_id, !!value) %>%
    as.data.frame()

  nObj <- summarized %>%
    distinct(group_id, nObject) %>%
    arrange(group_id) %>%
    pull(nObject)

  old_groups <- as.numeric(colnames(df)[-1])

  new_groups <- x$group_membership %>%
    arrange(match(old_groups, old_group)) %>%
    pull(new_group)

  # Color columns by group for easy viewing
  colors <- gen_color_palette(length(unique(new_groups)), output="colors")
  names(colors) <- unique(new_groups)

  colnames(df)[-1] <- sapply(2:ncol(df), function(idx){
    paste0("Group ", colnames(df)[idx], " (n=", nObj[[idx-1]], ")")
  })
  group_row <- c("New Group:", new_groups)
  df <- rbind(df, group_row)

  if(is_latex_output()){
    for (i in 2:ncol(df)){
      df[,i]<- cell_spec(df[,i], color="black", background=colors[[new_groups[[i-1]]]])
    }
    output <- kable(df, escape = F) %>%
      kable_styling(latex_options = "scale_down") %>%
      row_spec(nrow(df), color="black", bold=TRUE, background = "white", hline_after = TRUE) %>%
      knit_print()
  }else{
    output <- kable(df) %>% kable_styling()
    for (i in 2:ncol(df)){
      output <- column_spec(output, i, color="black", background=colors[[new_groups[[i-1]]]], border_left = TRUE, border_right = TRUE)
    }
    output <- row_spec(output, nrow(df), color="darkgray", bold=TRUE, background = "white", hline_after = TRUE)
  }
  knit_print(output)
}

#' @importFrom crayon make_style style
#' @export
print.cqc_cluster_panel <- function(x, ...){
  check_res <- attr(x, "check_res")
  value <- "marker"

  summarized <- check_res %>% summary

  df <- summarized %>%
    select(-nObject) %>%
    arrange(group_id) %>%
    spread(group_id, !!value) %>%
    as.data.frame()

  nObj <- summarized %>%
    distinct(group_id, nObject) %>%
    arrange(group_id) %>%
    pull(nObject)

  old_groups <- as.numeric(colnames(df)[-1])

  new_groups <- x$group_membership %>%
    arrange(match(old_groups, old_group)) %>%
    pull(new_group)

  # Color columns by group for easy viewing
  colors <- gen_color_palette(length(unique(new_groups)), output="styles")
  names(colors) <- unique(new_groups)

  channels <- df[,1]
  colorized <- data.frame(lapply(2:ncol(df), function(col){
    colors[[new_groups[[col-1]]]](df[,col])
  }))
  colorized <- cbind(df[,1], colorized)
  colnames(colorized) <- colnames(df)
  # Kable doesn't handle format decoration well. It bases header width on the
  # field width of the strings assessed with all formatting decoration (so super wide). So we need
  # to manually fix that
  rough <- kable(colorized, format = "rst")
  # Determining proper width of each column
  max.widths <- sapply(1:length(df), function(col){
    max(as.numeric(nchar(colnames(colorized)[col])), max(sapply(df[,col], nchar), na.rm = TRUE))
  })
  # Fix === rows
  headers <- lapply(max.widths, function(width){
    paste0(strrep("=", width), "  ")
  })
  headers <- do.call(paste0, headers)
  rough[[1]] <- rough[[3]] <- rough[[(length(rough))]] <- headers
  # Fix actual header row
  unpadded <- colnames(colorized)
  padded_top <- lapply(1:length(unpadded), function(idx){
    paste0(unpadded[[idx]], strrep(" ", max(0, max.widths[[idx]] - nchar(unpadded[[idx]]) + 2)))
  })
  padded_bottom <- lapply(1:length(new_groups), function(idx){
    paste0(new_groups[[idx]], strrep(" ", max(0, max.widths[[idx+1]] - nchar(new_groups[[idx]]) + 2)))
  })
  padded_bottom <- paste0("New Group:", strrep(" ", max(0, max.widths[[1]] - 8)), do.call(paste0, padded_bottom))

  rough[[2]] <- do.call(paste0, padded_top)

  buffer <- strrep(" ", (nchar(rough[[1]]) %/% 2) - 5)
  top <- paste0(buffer, "Old Groups", buffer)

  count_row <- lapply(1:length(nObj), function(idx){
    paste0(nObj[[idx]], strrep(" ", max(0, max.widths[[idx+1]] - nchar(nObj[[idx]]) + 2)))
  })
  count_row <- do.call(paste0, count_row)
  count_row <- paste0("n =", strrep(" ", max(0, max.widths[[1]] - 1)), count_row)

  rough <- c(top, rough, padded_bottom, headers, count_row)

  final <- structure(rough, class="knitr_kable", format="rst")
  print(final)
}

#' @export
knit_print.cqc_cluster <- function(x, ...){
  type <- gsub("cqc_cluster_", "", class(x)[[1]])
  type <- as.symbol(type)
  check_res <- as_tibble(attr(x, "check_res"))

  summarized <- check_res %>%
    inner_join(x$group_membership, by=c("group_id" = "old_group")) %>%
    select(-c(object)) %>%
    distinct() %>%
    group_by(group_id, new_group, nObject) %>%
    summarise(!!type := paste(!!type, collapse = ", ")) %>%
    arrange(group_id)

  new_groups <- x$group_membership$new_group

  colors <- gen_color_palette(length(unique(new_groups)), output="colors")
  names(colors) <- unique(new_groups)

  summarized <- kable(summarized) %>% kable_styling()

  # Should try to convert this for loop if possible, but row_spec makes it a little difficult
  for(i in 1:length(new_groups)){
    summarized <- row_spec(summarized, i, color = "black", background = colors[[new_groups[[i]]]])
  }
  knit_print(summarized)
}

#' @export
print.cqc_cluster <- function(x, ...){
  type <- gsub("cqc_cluster_", "", class(x)[[1]])
  check_res <- as_tibble(attr(x, "check_res"))
  summarized <- check_res %>%
    select(-c(object)) %>%
    distinct()

  type <- as.symbol(type)
  summarized <- group_by(summarized, group_id, nObject) %>%
    summarise(!!type := paste(!!type, collapse = ", ")) %>%
    arrange(group_id) %>%
    as.data.frame()

  new_groups <-  x$group_membership %>%
    arrange(old_group) %>%
    pull(new_group)

  summarized$new_group <- new_groups
  summarized <- summarized[,c(1,4,2,3)]

  # Color columns by group for easy viewing
  colors <- gen_color_palette(length(unique(new_groups)), output="styles")
  names(colors) <- unique(new_groups)

  colored_markers <- apply(summarized, 1, function(row){
    style(row[["marker"]], colors[[new_groups[[as.numeric(row[["group_id"]])]]]])
  })
  summarized$marker <- colored_markers

  final <- kable(summarized, format="rst")
  print(final)
}

object_type <- function(x){
  dat <- attr(x, "data")
  if(is(dat, "cqc_gs_list"))
    "GatingSet"
  else
    "FCS"
}

#' @export
print.cqc_check_summary <- function(x, collapse = TRUE, ...){
  if(collapse)
    x <- collapse_params(x)
  type <- object_type(x)
  colnames(x)[match("nObject", colnames(x))] = paste0("n", type)
  class(x) <- class(x)[-(1:3)]
  print(x)
}

#' @export
knit_print.cqc_check <- function(x, ...){
  knit_print(summary(x), ...)
}
#' Customized knit print for cqc_check_summary
#'
#'
#' @param x cqc_check_summary object returned by 'summary' call on `cqc_check`
#' @param collapse whether to collapse the same information within each group
#' @param ... not used
#' @importFrom dplyr ungroup everything
#' @export
knit_print.cqc_check_summary <- function(x, collapse = TRUE, ...) {
  n <- nrow(x)
  if (collapse) {
    x <- collapse_params(x)
  }
  type <- object_type(x)
  type <-  paste0("n", type)
  colnames(x)[match("nObject", colnames(x))] = type

  collaspse_idx <- match("group_id", colnames(x))
  if (is(x, "cqc_check_panel")) {
    collaspse_idx <- c(collaspse_idx, match(type, colnames(x)))
  }

  x <- x %>%
    arrange(desc(get(type))) %>%
    kable() %>%
    kable_styling("bordered", full_width = F, position = "left", font_size = 12)
  # browser()
  if (n > 0) {
    x <- x %>%
      collapse_rows(columns = collaspse_idx, "top") %>%
      row_spec(0, background = "#e5f5e0", color = "black")
  }

  knit_print(x)
}



#' @export
#' @importFrom colortable color_vctr
format.cqc_match_result_panel <- function(x, ...){
  class(x) <- class(x)[-(1:2)]
  res <- format(x, ...)
  #add color to ref
  res[[1]] <- color_vctr(res[[1]], text_color = "green")
  #find the ref col
  ref <- attr(x, "ref")
  ref_col <- paste0("group ", ref)
  ref_idx <- grep(ref_col, colnames(res), fixed = TRUE)
  res[[ref_idx]] <- color_vctr(res[[ref_idx]], text_color = "green")
  colnames(res)[ref_idx] <- "Ref group"

  res
}

#' @export
#' @importFrom tidyr spread
format.cqc_check_panel <- function(x, color_ref = FALSE, ...){
  # x <- summary(x)
  anchor <- attr(x, "by")
  # browser()
  if(anchor == "channel")
    value <- "marker"
  else
    value <- "channel"
  # #long to wide
   x %>% summary %>%
    mutate(group_id := paste("group", group_id), nObject := paste0("(n=", nObject, ")")) %>%
    unite(grp, group_id, nObject, sep = "") %>% #merge grp cols
    spread(grp, !!value) %>%
    filter(get(anchor) !="") %>% #rm the empty row that was caused by samples that have entire empty markers
    `class<-`(value = class(x)[-(1:2)])

}

#' @export
print.cqc_match_result_panel <- function(x, ...){
  format(x, ...) %>%
    print
}

#' @export
print.cqc_check_panel <- function(x, ...){
  format(x, ...) %>%
    print
}

#' @export
knit_print.cqc_match_result_panel <- function(x, ...){

  x <- kable_cqc_check_panel(x, ...)
  ref_idx <- attr(x, "ref")
  x %>%
    column_spec(c(1, ref_idx), color = "green") %>% knit_print

}

#' @export
knit_print.cqc_check_panel <- function(x, ...){
  kable_cqc_check_panel(x, ...)%>%
    knit_print

}
#' @importFrom dplyr mutate_all
kable_cqc_check_panel <- function(x, ...){

  attr <- attributes(x)

  x <- mutate_all(x, htmlEscape)
  attributes(x) <- attr
  x <- format(x, ...)

  x[is.na(x)] <- "<font color='red'>N/A</font>"
  ref_idx <- grep("Ref", colnames(x), fixed = TRUE)

  # x[x=="FSC-Height"] <- "<font color='red'>N/A</font>"
  x = x%>%
    kable(escape = F) %>%
    kable_styling(c("bordered", "condensed"), full_width = F, position = "left", font_size = 12)
  attr(x, "ref") <- ref_idx
  x
}
gen_color_palette <- function(n, output=c("styles", "colors")){
  output <- match.arg(output, c("styles", "colors"))
  colors <- rainbow(n, s = 0.3)
  if(output == "colors")
    colors
  else
    lapply(colors, function(color){
      make_style(color)
    })
}
