

#' Extract channel/marker info from a cytoframe object
#'
#' @param cf cytoframe
#' @param skip_na whether to skip the entries with marker unset (i.e. NA)
#'
#' @return a tibble with two columns: "channel" and 'marker'
#' @importFrom dplyr select rename
#' @importFrom flowCore parameters
#' @examples 
#' fcs_path <- system.file("extdata", "GvHD_QC", "s5a01.fcs", package = "cytoqc")
#' cf <- load_cytoframe_from_fcs(fcs_path)
#' cf_get_panel(cf)
#' @export
cf_get_panel <- function(cf, skip_na = FALSE) {
  res <- pData(parameters(cf)) %>%
    as_tibble() %>%
    select(c("name", "desc")) %>%
    rename(channel = name, marker = desc)

  if(skip_na)
    res <- filter(res, is.na(marker) == FALSE)
  #remove AsIs attribute
  class(res$channel) <- NULL
  class(res$marker) <- NULL

  res
}

#' @rdname cqc_check
#' @export
cqc_check_channel <- function(x, ...){
  cqc_check(x, type = "channel", ...)
}
#' @rdname cqc_check
#' @export
cqc_check_marker <- function(x, ...){
  cqc_check(x, type = "marker", ...)
}

#' @param  by the column used as the anchor when 'x' is the check result of panel, it can be either 'channel' or 'marker'
#' @rdname cqc_check
#' @export
cqc_check_panel <- function(x, by = "channel", ...){
  cqc_check(x, type = "panel", by = by, ...)
}
#' @rdname cqc_check
#' @export
cqc_check_keyword <- function(x, ...){
  cqc_check(x, type = "keyword", ...)
}
#' @rdname cqc_check
#' @export
cqc_check_gate <- function(x, ...){
  cqc_check(x, type = "gate", ...)
}

#' Perform a QC check on flow data.
#'
#' This is the first step of the entire qc workflow.
#' It extracts meta information(specified by 'type' argument) from the raw data
#' and compare/group them across samples.
#' This provides a sample-wise data table for the further summary report.
#'
#' @return a tibble with 4 columns: object, qc type (e.g. channel), group_id and nobject (i.e. group count)
#' @param x \code{\link{cqc_cf_list}}, \code{\link{cqc_gs}}, or \code{\link{cqc_gs_list}} object
#' @param ... additional arguments.
#' 
#'  type -- specify the qc type, can be "channel", "marker" or "panel"
#'
#'  delimiter -- a special character used to separate channel and marker
#'
#'  keys -- The vector to supply the keys to be grouped on. default is NULL, which is extracted automatically from the flow data
#'
#'
#' @examples
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' 
#' # You may directly call the method for the parameter you would like to check
#' keyword_groups <- cqc_check_keyword(qc_cf_list)
#' keyword_groups
#' 
#' # Or use the type argument
#' channel_groups <- cqc_check(qc_cf_list, type = "channel")
#' channel_groups
#' 
#' panel_groups <- cqc_check(qc_cf_list, type = "panel", by = "marker")
#' panel_groups
#' 
#' @export
cqc_check <- function(x, ...) UseMethod("cqc_check")

#' @export
cqc_check.cqc_gs_list <- function(x, type, delimiter = "|", ...) {
  type = match.arg(type, c("channel", "marker", "panel", "keyword", "gate"))
  #extract the first cf from each gs
  cflist <- sapply(x, function(gs) get_cytoframe_from_cs(gs_pop_get_data(gs), 1)) # TODO:qc within gs to ensure all data are consistent
  cflist <- cqc_cf_list(cflist)
  if (type == "gate") {
    sep <- paste0(delimiter, delimiter)
    keys <- sapply(x, function(gs) {
      key <- gs_get_pop_paths(gs, path = "auto")#extract gate paths
      #order them alphabetically and collapse them as a string
      key %>%
        sort() %>%
        paste(collapse = sep)
    })
  } else {
    keys <- NULL
  }
  #dispatch to cqc_check.cqc_cf_list
  res <- cqc_check(cflist, type, keys, delimiter, ...)
  # class(res) <- c("cqc_check_gs", class(res))

  attr(res, "data") <- x
  res
}

#' @export
#' @importFrom dplyr filter arrange pull mutate group_indices distinct count add_count
#' @importFrom tidyr separate separate_rows
cqc_check.cqc_cf_list <- function(x, type, keys = NULL, delimiter = "|", by = "channel", ...) {
  by <- match.arg(by, c("channel", "marker"))
  types <-  c("channel", "marker", "panel", "keyword")
  if(!is.null(keys))#passed down from cqc_check.cqc_gs_list
    types <- c(types, "gate")
  type = match.arg(type, types)
  sep <- paste0(delimiter, delimiter) # double delimiter for sep params and single delimiter for sep channel and marker
  #if keys are not supplied then extract them from data according to the type
   if (is.null(keys)) {
    keys <- sapply(x, function(cf) {
      if (type == "keyword") {
        key <- names(keyword(cf, compact = TRUE))
      } else {
        key <- cf_get_panel(cf) # %>% arrange(channel)
        if (type != "channel") {
          key <- filter(key, is.na(marker) == FALSE)
        }
        if (type == "panel") {
          key <- unite(key, panel, channel, marker, sep = delimiter)
        }
        key <- key[[type]]
      }

      key %>%
        sort() %>%
        paste(collapse = sep)
    })
  }
  #convert to itemized(one entry per object&key) tbl
  res <- tibble(object = names(keys), key = keys)
  #generate group id based on the key
  gid <- res %>%
    group_by(key) %>%
    group_indices()
  res <- res %>%
    mutate(group_id = gid) %>%
    add_count(group_id, key) %>%
    rename(nObject = n) %>%
    separate_rows(key, sep = paste0("\\Q", sep, "\\E"))#split the collapsed key string into separate rows(one key per row)
  if (type == "panel") {#for panel check, need to split chnl and marker into separate columns
    res <- separate(res, key, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
    attr(res, "by") <- by

  } else {
    res <- rename(res, !!type := key)
  }


  class(res) <- c("cqc_check", class(res))
  class(res) <- c(paste0("cqc_check_", type), class(res))
  attr(res, "data") <- x
  res
}

#' @export
cqc_check.cqc_gs <- function(x, type, keys = NULL, delimiter = "|", ...) {
  type = match.arg(type, c("channel", "marker", "panel", "keyword", "gate"))
  # Get list of cytoframes for dispatching to cqc_cflist methdods
  cflist <- lapply(1:length(x), function(idx) {gh_pop_get_data(x[[idx]], returnType="cytoframe")})
  names(cflist) <- names(x)
  cflist <- cqc_cf_list(cflist)

  if (type == "gate") {
    sep <- paste0(delimiter, delimiter)
    keys <- sapply(x, function(gh) {
      key <- gh_get_pop_paths(gh, path = "auto")#extract gate paths
      #order them alphabetically and collapse them as a string
      key %>%
        sort() %>%
        paste(collapse = sep)
    })
  } else {
    keys <- NULL
  }

  #dispatch to cqc_check.cqc_cf_list
  res <- cqc_check(cflist, type, keys, delimiter, ...)

  attr(res, "data") <- x
  res
}

#' Provide the summary view of the cqc_check report.
#'
#' It summarize the itemized check report into group overview.
#' It is automatically triggered by print method of `cqc_check` result.
#'
#' @param object qc table returned by \code{\link{cqc_check}}
#' @param ... Additional arguments not for the user. Ignore.
#' @examples
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' 
#' channel_groups <- cqc_check(qc_cf_list, type = "channel")
#' summary(channel_groups)
#' 
#' @export
summary.cqc_check <- function(object, ...) {
  res <- object %>%
    select(-c(object)) %>%
    distinct()
  class(res) <- c("cqc_check_summary", class(res))
  attr(res, "data") <- attr(object, "data")
  res
}

#' @export
diff.cqc_check_gate <- function(x, ...) {
  diff.cqc_check(x, c("gate"))
}
#' @export
diff.cqc_check_keyword <- function(x, ...) {
  diff.cqc_check(x, c("keyword"))
}

#' @export
diff.cqc_check_channel <- function(x, ...) {
  diff.cqc_check(x, c("channel"))
}

#' @export
diff.cqc_check_marker <- function(x, ...) {
  diff.cqc_check(x, c("marker"))
}
#' @export
diff.cqc_check_panel <- function(x, ...) {
  diff.cqc_check(x, c("channel", "marker"))
}

#' Display the differences among QC groups
#'
#' @param x qc table returned by \code{\link{cqc_check}}
#' @param vars variable to split by. Determined automatically.
#' @param ... Additional arguments not for the user. Ignore.
#' @examples
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' 
#' channel_groups <- cqc_check(qc_cf_list, type = "channel")
#' diff(channel_groups)
#' @importFrom dplyr group_split inner_join anti_join
#' @importFrom purrr reduce map map_dfr
#' @export
diff.cqc_check <- function(x, vars, ...) {
  grps <- x %>% summary %>%
    group_split(group_id)
  commons <- grps %>% reduce(inner_join, by = vars)
  grps %>%
    map_dfr(anti_join, y = commons, by = vars) %>%
    `class<-`(value =c("cqc_check_summary", class(x))) %>%
    `attr<-`("data", value = attr(x, "data"))

}

#' Split the result of \code{cqc_check} into groups
#'
#' It is used to split samples into separate groups when they can't be reconciled into the same group.
#'
#' @importFrom purrr walk
#' @param x cqc_check object
#' @param f,drop,... not used
#' @examples 
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' channel_groups <- cqc_check(qc_cf_list, type = "channel")
#' # Show the different groups that will result from the split
#' summary(channel_groups)
#' split_groups <- split(channel_groups)
#' split_groups
#' @export
split.cqc_check <- function(x, f, drop = FALSE, ...) {
  cqc_data <- attr(x, "data")
  data_type <- class(cqc_data)
  vec <- x %>%
    select(c(object, group_id)) %>%
    distinct() %>%
    pull(group_id)
  split(cqc_data, vec) %>% map(function(i) {
    class(i) <- c(data_type, class(i))
    i
  })
}

#' visualize the tree structure difference among the GatingSets
#'
#' @param groups \code{cqc_check_gate} grouping resulte from \code{cqc_check}.
#' @examples
#' gs_paths <- list.files(system.file("extdata", "gslist_manual_QC", package = "cytoqc"), full.names = TRUE)
#' gs1 <- load_gs(gs_paths[[1]])
#' gs2 <- load_gs(gs_paths[[2]])
#' qc_gslist <- cqc_gs_list(list(gs1, gs2))
#' groups <- cqc_check(qc_gslist, type="gate")
#' plot_diff(groups)
#' 
#' @export
#' @import Rgraphviz graph
plot_diff <- function(groups) {
  grps <- split(groups)
  df <- diff((groups))

  nGrps <- length(grps)
  # par(mfrow = c(1,nGrps))

  for (grpid in names(grps))
  {
    message("processing group: ", grpid)
    gs <- grps[[grpid]][[1]]
    if (!is(gs, "GatingSet")) {
      stop("Expecting a GatingSet!")
    }
    gh <- gs[[1]]

    g <- flowWorkspace:::.getGraph(gh)
    nodes.uncommon <- filter(df, group_id == grpid) %>% pull(gate)
    nodeDataDefaults(g, "common") <- FALSE

    message("hide/tag the common nodes ...")
    node.path2 <- gs_get_pop_paths(gh, showHidden = TRUE, path = 3)
    thisNodeSet <- gs_get_pop_paths(gs, path = "auto")
    commonNodes <- setdiff(thisNodeSet, nodes.uncommon)

    for (node in commonNodes)
    {
      nodeID <- flowWorkspace:::.getNodeInd(gh, node) - 1

      node.label <- paste0("N_", nodeID)
      nodeData(g, node.label, attr = "common") <- TRUE
      children <- gs_pop_get_children(gs[[1]], node, path = "auto")
      # keep the direct parent of uncommon node
      if (!any(children %in% nodes.uncommon) || node == "root") {
        nodeData(g, node.label, attr = "hidden") <- "1"
      } else {
        nodeData(g, node.label, attr = "label") <- node.path2[nodeID + 1]
      }
    }


    # filter out hidden node
    nodes.hidden <- nodeData(g, attr = "hidden")
    for (i in 1:length(nodes.hidden))
    {
      if (as.logical(as.integer(nodes.hidden[[i]]))) {
        nodeID <- names(nodes.hidden[i])

        g <- removeNode(nodeID, g)
      }
    }
    message("Rendering the substree ...")
    graphAttr <- list(rankdir = "LR", page = c(8.5, 11))
    nAttrs <- list()
    nodes <- nodes(g)
    if (length(nodes) == 0) {
      g <- addNode("Root", g)
      g <- addNode("Common nodes", g)
      g <- addEdge("Root", "Common nodes", g)
      nodes <- nodes(g)
      nLabel <- nodes
      names(nLabel) <- nodes
      plot(g,
        nodeAttrs = list(label = nLabel),
        attrs = list(
          graph = graphAttr,
          node = list(
            fixedsize = FALSE,
            fillcolor = "gray",
            shape = "rectangle"
          )
          # ,edge = list(col = "transparent")
        ),
        main = paste0("Group ", grpid)
      )
      # plot.new()
      next
    }

    common.flags <- nodeData(g, attr = "common")
    nAttrs$label <- unlist(nodeData(g, attr = "label"))
    nAttrs$fillcolor <- sapply(
      common.flags,
      function(iscommon) {
        ifelse(iscommon, "gray", "red")
      }
    )
    nAttrs$lty <- sapply(
      common.flags,
      function(iscommon) {
        ifelse(!iscommon, "dotted", "solid")
      }
    )

    # pass plot parameters to node attributes (some of parameters won't work via passing to layoutGraph directly)
    nAttrs[["fixedsize"]] <- sapply(common.flags, function(i) FALSE)
    nAttrs[["shape"]] <- sapply(
      common.flags,
      function(iscommon) {
        ifelse(iscommon, "rectangle", "ellipse")
      }
    )
    # params <- list(...)
    # for(pname in names(params))
    #   nAttrs[[pname]] <- sapply(nodes, function(i)params[[pname]])
    nodeRenderInfo(g) <- nAttrs

    plot(g,
      nodeAttrs = nAttrs, attrs = list(graph = graphAttr),
      main = paste0("Group ", grpid)
    )
  }
}

#' Remove outlier groups from the result of 'cqc_check'.
#'
#' @param groups the object returned by 'cqc_checks'
#' @param id the group id to be dropped from the dataset
#' @examples 
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' channel_groups <- cqc_check(qc_cf_list, type = "channel")
#' summary(channel_groups)
#' channel_groups <- cqc_drop_groups(channel_groups, 2)
#' channel_groups
#' @export
cqc_drop_groups <- function(groups, id) {
  cqc_data <- attr(groups, "data")
  data_type <- class(cqc_data)
  torm <- filter(groups, group_id == id) %>%
    pull(object) %>%
    unique()
  cqc_data <- cqc_data[-match(torm, names(cqc_data))]
  class(cqc_data) <- data_type
  groups <- filter(groups, group_id != id)
  attr(groups, "data") <- cqc_data
  groups
}

#' Extract the data from the result of a \code{cqc_check} call.
#'
#' @param groups the object returned by \code{\link{cqc_checks}}
#' @param id the group id to be selected from the dataset, default is NULL, meaning all data
#' @examples 
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' channel_groups <- cqc_check(qc_cf_list, type = "channel")
#' summary(channel_groups)
#' group_3_cf_list <- cqc_get_data(channel_groups, 3)
#' # A list of cytoframes
#' group_3_cf_list
#' @export
cqc_get_data <- function(groups, id = NULL) {
  cqc_data <- attr(groups, "data")
  if (!is.null(id)) {
    sel <- filter(groups, group_id == id) %>%
      pull(object) %>%
      unique()
    cqc_data <- cqc_data[sel]
  }
  cqc_data
}

