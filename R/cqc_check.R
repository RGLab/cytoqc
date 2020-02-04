

#' Extract channel/marker info from a cytoframe object
#'
#' @param cf cytoframe
#' @param skip_na whether to skip the entries with marker unset (i.e. NA)
#'
#' @return a tibble with two columns: "channel" and 'marker'
#' @importFrom tibble as.tibble
#' @importFrom dplyr select rename
#' @importFrom flowCore parameters
#' @export
cf_get_panel <- function(cf, skip_na = FALSE) {
  res <- pData(parameters(cf)) %>%
    as.tibble() %>%
    select(c("name", "desc")) %>%
    rename(channel = name, marker = desc)

  if(skip_na)
    res <- filter(res, is.na(marker) == FALSE)

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
#' @rdname cqc_check
#' @export
cqc_check_panel <- function(x, ...){
  cqc_check(x, type = "panel", ...)
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
#' @param x cqc_cf_list
#' @param ... additional arguments.
#'  type specify the qc type, can be "channel", "marker" or "panel"
#'  delimiter a special character used to separate channel and marker
#'  keys The vector to supply the keys to be grouped on. default is NULL, which is extracted automatically from the flow data
#' @examples
#' \dontrun{
#' groups <- cqc_check(cqc_cf_list, "channel")
#' }
#' @export
cqc_check <- function(x, ...) UseMethod("cqc_check")

#' @export
cqc_check.cqc_gs_list <- function(x, type, delimiter = "|", ...) {
  type = match.arg(type, c("channel", "marker", "panel", "keyword", "gate"))
  cflist <- sapply(x, function(gs) get_cytoframe_from_cs(gs_pop_get_data(gs), 1)) # TODO:qc within gs to ensure all data are consistent
  cflist <- cqc_cf_list(cflist)
  if (type == "gate") {
    sep <- paste0(delimiter, delimiter)
    keys <- sapply(x, function(gs) {
      key <- gs_get_pop_paths(gs, path = "auto")

      key %>%
        sort() %>%
        paste(collapse = sep)
    })
  } else {
    keys <- NULL
  }
  res <- cqc_check(cflist, type, keys, delimiter, ...)
  # class(res) <- c("cqc_check_gs", class(res))

  attr(res, "data") <- x
  res
}

#' @export
#' @importFrom dplyr filter arrange pull mutate group_indices distinct count add_count
#' @importFrom tidyr separate separate_rows
cqc_check.cqc_cf_list <- function(x, type, keys = NULL, delimiter = "|", ...) {
  types <-  c("channel", "marker", "panel", "keyword")
  if(!is.null(keys))#passed down from cqc_check.cqc_gs_list
    types <- c(types, "gate")
  type = match.arg(type, types)
  sep <- paste0(delimiter, delimiter) # double delimiter for sep params and single delimiter for sep channel and marker
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

  res <- tibble(object = names(keys), key = keys)
  gid <- group_indices(res, key)
  res <- res %>%
    mutate(group_id = gid) %>%
    add_count(group_id, key) %>%
    rename(nObject = n) %>%
    separate_rows(key, sep = paste0("\\Q", sep, "\\E"))
  if (type == "panel") {
    res <- separate(res, key, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
  } else {
    res <- rename(res, !!type := key)
  }


  #
  # res <- strsplit(res, split= sep, fixed = "TRUE")[[1]]
  # res <- tibble(reference = res)
  # if(type == "marker")
  #   res <- separate(res, channel, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
  class(res) <- c("cqc_check", class(res))
  class(res) <- c(paste0("cqc_check_", type), class(res))
  attr(res, "data") <- x
  res
}

#' Provide the summary view of the cqc_check report.
#'
#' It summarise the sample-wise qc report into group overview report for reference selection and further QC action.
#'
#' @param object qc table returned by 'cqc_check'
#' @param ... Additional arguments not for the user. Ignore.
#' @examples
#' \dontrun{
#' su <- summary(groups)
#' }
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
#' @param x the grouped summary report generated by 'summary' call on the 'cqc_check' results
#' @param vars variable to split by. Determined automatically.
#' @param ... Additional arguments not for the user. Ignore.
#' @examples
#' \dontrun{
#' su <- summary(groups)
#' diff(su)
#' }
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

#' Split the result of 'cqc_check' into groups
#'
#' It is used to split samples into separate groups when they can't be reconciled into the sampe group.
#'
#' @importFrom purrr walk
#' @param x cqc_check object
#' @param f,drop,... not used
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

#' visualize the tree structure differnece among the GatingSets
#'
#' @param groups \code{cqc_check_gate} grouping resulte from \code{cqc_check}.
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

#' Extract the data from the result of a 'cqc_checks' call.
#'
#' @param groups the object returned by 'cqc_checks'
#' @param id the group id to be selected from the dataset, default is NULL, meaning all data
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
#' Set a reference sample for further QC.
#'
#' This is the first step prior to finding a QC solution for a data set.
#'
#' @param x cqc report generated by 'cqc_check'
#' @param ref specifies the reference, which can be either an integer group id or a characte vector giving the actual values of the reference
#' @return the original cqc report with the reference info attached
cqc_set_reference <- function(x, ref) {
  .Defunct("cqc_match")
  attr(x, "reference") <- ref
  x
}
