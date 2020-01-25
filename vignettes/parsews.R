library(flowWorkspace)
library(CytoML)

path <- "~/remote/fh/fast/gottardo_r/mike_working/lyoplate_out/"

#' flowJo workspace path
ws_path <- file.path(path, "XML")
#' raw FCS file path
fcs_path <- file.path(path, "FCS")
#' path that sores the phenoData as excel file
pd_path <- file.path(path, "pData")
out_path <- file.path(path, "parsed")


#centers <- list.dirs(ws_path, recursive = FALSE, full.names = FALSE)
centers <- c('BIIR','CIMR','Miami','NHLBI','Stanford','UCLA','Yale')

#' We parse each center's workspace and then extract its flowSet.
message("Parsing Lyoplate workspaces...")
for(center in centers) {
  message("Center: ", center)
  this_out <- file.path(out_path, center)
  if(!dir.exists(this_out))
    dir.create(this_out)
  this_ws_path <- file.path(ws_path, center, paste0("CA c3 v2 ", center, ".xml"))
  ws <- open_flowjo_xml(this_ws_path)
  groupNames <- levels(fj_ws_get_sample_groups(ws)[,"groupName"])
  #NOTE: here we use the fussy match due to the short name used in auto gating
  #ideally, we want to do the strict full string match avoid picking up the wrong panel
  #in case the panel name are similar. But in this particular panel set, it is safe to do so.
  for(panel in c("tcell", "bcell", "DC", "treg"))
  {
    this_gs_path <- file.path(this_out, panel)
    if(!dir.exists(this_gs_path))
    {
      message("parsing: ", panel)

      groupID <- grep(panel, groupNames, ignore.case = T)
      #NOTE: in theory, we shouldn't need to do this if center names were consistent.
      this_fcs_path <- file.path(fcs_path, center)
      this_fcs_path <- sub("BIIR", "Baylor", this_fcs_path)
      gs <- flowjo_to_gatingset(ws, name = groupNames[groupID], path = this_fcs_path, extend_val = -10)
      save_gs(gs, this_gs_path)
    }


  }

}





