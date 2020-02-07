
# QC for tcell panel

``` r
library(flowCore)
library(flowWorkspace)
library(cytoqc)
# library(printr)
# library(DT)
```

``` r
path <- "~/remote/fh/fast/gottardo_r/mike_working/lyoplate_out/parsed"
centers <- c('BIIR','CIMR','Miami','NHLBI','Stanford','UCLA','Yale')
```

## Load gs

``` r
panel <- "tcell"
gslist <- sapply(centers, function(center) {
  message("Center: ", center)
  gs <- load_gs(file.path(path, center, panel))
})

cqc_data <- cqc_gs_list(gslist)
```

## QC Check gates

``` r
#group by gates
groups <- cqc_check(cqc_data, "gate")
groups
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

gate

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

4- 8+, 4- 8+/38- DR-, 4- 8+/38- DR+, 4- 8+/38+ DR-, 4- 8+/38+ DR+, 4-
8+/CCR7- 45RA-, 4- 8+/CCR7- 45RA+, 4- 8+/CCR7+ 45RA-, 4- 8+/CCR7+ 45RA+,
4+ 8-, 4+ 8-/38- DR-, 4+ 8-/38- DR+, 4+ 8-/38+ DR-, 4+ 8-/38+ DR+, 4+
8-/CCR7- 45RA-, 4+ 8-/CCR7- 45RA+, 4+ 8-/CCR7+ 45RA-, 4+ 8-/CCR7+ 45RA+,
CD3, DNT, DPT, LYM, not dead, root, singlets

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

4- 8+, 4- 8+/38- DR-, 4- 8+/38- DR+, 4- 8+/38+ DR-, 4- 8+/38+ DR+, 4-
8+/CCR7- 45RA-, 4- 8+/CCR7- 45RA+, 4- 8+/CCR7+ 45RA-, 4- 8+/CCR7+ 45RA+,
4+ 8-, 4+ 8-/38- DR-, 4+ 8-/38- DR+, 4+ 8-/38+ DR-, 4+ 8-/38+ DR+, 4+
8-/CCR7- 45RA-, 4+ 8-/CCR7- 45RA+, 4+ 8-/CCR7+ 45RA-, 4+ 8-/CCR7+ 45RA+,
CD3, DNT, DPT, LYM, not dead,
root

</td>

</tr>

</tbody>

</table>

``` r
diff(groups)
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

gate

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

singlets

</td>

</tr>

</tbody>

</table>

``` r
#vis the difference
# plot_diff(groups)

# match reference
match_result <- cqc_match(groups, ref = 1)
# match_result

cqc_fix(match_result)

cqc_check(cqc_data, "gate")
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

gate

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

4- 8+, 4- 8+/38- DR-, 4- 8+/38- DR+, 4- 8+/38+ DR-, 4- 8+/38+ DR+, 4-
8+/CCR7- 45RA-, 4- 8+/CCR7- 45RA+, 4- 8+/CCR7+ 45RA-, 4- 8+/CCR7+ 45RA+,
4+ 8-, 4+ 8-/38- DR-, 4+ 8-/38- DR+, 4+ 8-/38+ DR-, 4+ 8-/38+ DR+, 4+
8-/CCR7- 45RA-, 4+ 8-/CCR7- 45RA+, 4+ 8-/CCR7+ 45RA-, 4+ 8-/CCR7+ 45RA+,
CD3, DNT, DPT, LYM, not dead, root

</td>

</tr>

</tbody>

</table>

## QC check for channel

``` r
groups <- cqc_check(cqc_data, "channel")
groups
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

channel

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<Alexa Fluor 488-A\>, \<AmCyan-A\>, \<APC-A\>, \<APC-Cy7-A\>, \<Pacific
Blue-A\>, \<PE YG-A\>, \<PE-Cy7 YG-A\>, \<PerCP-Cy5-5-A\>, FSC-A, SSC-A,
Time

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<Am Cyan-A\>, \<APC-A\>, \<APC-Cy7-A\>, \<FITC-A\>, \<Pacific Blue-A\>,
\<PE-A\>, \<PE-Cy7-A\>, \<PerCP-Cy5-5-A\>, FSC-A, FSC-H, FSC-W, SSC-A,
SSC-H, SSC-W, Time

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<AmCyan-A\>, \<APC-A\>, \<APC-Cy7-A\>, \<FITC-A\>, \<Pacific Blue-A\>,
\<PE Cy7 YG-A\>, \<PE-A\>, \<PerCP-Cy5-5-A\>, FSC-A, FSC-W, SSC-A,
SSC-W, Time

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<APC-A\>, \<APC-Cy7-A\>, \<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>,
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>, FSC-A, SSC-A, Time

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<APC-A\>, \<APC-H7-A\>, \<BD Horizon V450-A\>, \<BD Horizon V500-A\>,
\<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>, \<PerCP-Cy5-5-A\>, FSC-A, FSC-H,
FSC-W, SSC-A, SSC-H, SSC-W, Time

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<APC-A\>, \<APC-H7-A\>, \<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>,
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>, FSC-A, FSC-H, SSC-A, SSC-H,
Time

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<B515-A\>, \<B710-A\>, \<G560-A\>, \<G780-A\>, \<R660-A\>, \<R780-A\>,
\<V450-A\>, \<V545-A\>, FSC-A, FSC-H, FSC-W, SSC-A, SSC-H, SSC-W, Time

</td>

</tr>

</tbody>

</table>

## Channels are very different across centers so move on to check marker

``` r
groups <- cqc_check(cqc_data, "marker")
groups
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

marker

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CD197, CD3, CD38, CD4, CD45RA, CD8, HLA-DR, Live Green

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CCR7 PE, CD3 V450, CD38 APC, CD4 PerCP-Cy55, CD45RA PE-Cy7, CD8 APC-H7,
HLA-DR V500, Live Dead FITC

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CCR7, CD3, CD38, CD4, CD45RA, CD8, HLA DR, Live/Dead

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CCR7, CD3, CD38, CD4, CD45RA, CD8, HLA-DR, LIVE

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CCR7, CD3, CD38, CD4, CD45RA, CD8, HLADR, LIVE\_GREEN

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CD197, CD3, CD38, CD4, CD45RA, CD8, HLA-DR, LIVE DEAD

</td>

</tr>

</tbody>

</table>

## Markers are more standardized and go ahead to further clean it

``` r
res <- cqc_match(groups, ref = 3) 
res
```

## Re-match by relaxing the string matching threshold

``` r
res <- cqc_match(groups, ref = 3, max.distance = 0.6)
res
```

## Manually match the individual items that are still not matched

``` r
res <- cqc_update_match(res, group = 1, map = c("CD4 PerCP-Cy55" = "CD4"
                                              , "CD8 APC-H7" = "CD8"
                                              , "Live Dead FITC" = "LIVE")
                        )
res
```

``` r
cqc_fix(res)
```

## update checks

``` r
groups <- cqc_check(cqc_data, "marker")
groups
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

marker

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CCR7, CD3, CD38, CD4, CD45RA, CD8, HLA-DR, LIVE

</td>

</tr>

</tbody>

</table>

## check pannel

``` r
res <- cqc_check(cqc_data, "panel")
res
```

    ## # A tibble: 25 x 8
    ##    channel `group 1(n=1)` `group 2(n=1)` `group 3(n=1)` `group 4(n=1)`
    ##    <chr>   <chr>          <chr>          <chr>          <chr>         
    ##  1 <Alexa… LIVE           <NA>           <NA>           <NA>          
    ##  2 <Am Cy… <NA>           HLA-DR         <NA>           <NA>          
    ##  3 <AmCya… HLA-DR         <NA>           HLA-DR         <NA>          
    ##  4 <APC-A> CD38           CD38           CD38           CD38          
    ##  5 <APC-C… CD8            CD8            CD8            CD8           
    ##  6 <APC-H… <NA>           <NA>           <NA>           <NA>          
    ##  7 <B515-… <NA>           <NA>           <NA>           <NA>          
    ##  8 <B710-… <NA>           <NA>           <NA>           <NA>          
    ##  9 <BD Ho… <NA>           <NA>           <NA>           <NA>          
    ## 10 <BD Ho… <NA>           <NA>           <NA>           <NA>          
    ## # … with 15 more rows, and 3 more variables: `group 5(n=1)` <chr>, `group
    ## #   6(n=1)` <chr>, `group 7(n=1)` <chr>

### Spread markers

``` r
format(res, anchor = "marker")
```

    ## # A tibble: 8 x 8
    ##   marker `group 1(n=1)` `group 2(n=1)` `group 3(n=1)` `group 4(n=1)`
    ##   <chr>  <chr>          <chr>          <chr>          <chr>         
    ## 1 CCR7   <PE YG-A>      <PE-A>         <PE-A>         <PE-A>        
    ## 2 CD3    <Pacific Blue… <Pacific Blue… <Pacific Blue… <V450-A>      
    ## 3 CD38   <APC-A>        <APC-A>        <APC-A>        <APC-A>       
    ## 4 CD4    <PerCP-Cy5-5-… <PerCP-Cy5-5-… <PerCP-Cy5-5-… <PerCP-Cy5-5-…
    ## 5 CD45RA <PE-Cy7 YG-A>  <PE-Cy7-A>     <PE Cy7 YG-A>  <PE-Cy7-A>    
    ## 6 CD8    <APC-Cy7-A>    <APC-Cy7-A>    <APC-Cy7-A>    <APC-Cy7-A>   
    ## 7 HLA-DR <AmCyan-A>     <Am Cyan-A>    <AmCyan-A>     <V500-A>      
    ## 8 LIVE   <Alexa Fluor … <FITC-A>       <FITC-A>       <FITC-A>      
    ## # … with 3 more variables: `group 5(n=1)` <chr>, `group 6(n=1)` <chr>, `group
    ## #   7(n=1)` <chr>

## Use the marker as reference to standardize the channels across centers

``` r
cf <- gs_cyto_data(cqc_data[[1]])[[1]]
panel <- cf_get_panel(cf, skip_na = TRUE)
panel
```

    ## # A tibble: 8 x 2
    ##   channel         marker  
    ##   <I<chr>>        <I<chr>>
    ## 1 <APC-A>         CD38    
    ## 2 <APC-H7-A>      CD8     
    ## 3 <FITC-A>        LIVE    
    ## 4 <PerCP-Cy5-5-A> CD4     
    ## 5 <V450-A>        CD3     
    ## 6 <V500-A>        HLA-DR  
    ## 7 <PE-A>          CCR7    
    ## 8 <PE-Cy7-A>      CD45RA

``` r
cqc_set_panel(cqc_data, panel, ref.col = "marker")
groups <- cqc_check(cqc_data, "panel")
groups
```

<!--html_preserve-->

<div id="htmlwidget-d50c9c488f932329f67b" class="datatables html-widget" style="width:960px;height:500px;">

</div>

<script type="application/json" data-for="htmlwidget-d50c9c488f932329f67b">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8"],["&lt;APC-A&gt;","&lt;APC-H7-A&gt;","&lt;FITC-A&gt;","&lt;PE-A&gt;","&lt;PE-Cy7-A&gt;","&lt;PerCP-Cy5-5-A&gt;","&lt;V450-A&gt;","&lt;V500-A&gt;"],["CD38","CD8","LIVE","CCR7","CD45RA","CD4","CD3","HLA-DR"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>channel<\/th>\n      <th>group 1(n=7)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}]}},"evals":[],"jsHooks":[]}</script>

<!--/html_preserve-->

## Refresh QC report

``` r
groups <- cqc_check(cqc_data, "channel")
groups
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

channel

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

\<APC-A\>, \<APC-H7-A\>, \<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>,
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>, FSC-A, FSC-H, FSC-W, SSC-A,
SSC-H, SSC-W, Time

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

\<APC-A\>, \<APC-H7-A\>, \<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>,
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>, FSC-A, SSC-A, Time

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<APC-A\>, \<APC-H7-A\>, \<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>,
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>, FSC-A, FSC-H, SSC-A, SSC-H,
Time

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

\<APC-A\>, \<APC-H7-A\>, \<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>,
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>, FSC-A, FSC-W, SSC-A, SSC-W,
Time

</td>

</tr>

</tbody>

</table>

``` r
diff(groups)
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nGatingSet

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

channel

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

FSC-H, FSC-W, SSC-H, SSC-W

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

FSC-H, SSC-H

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

FSC-W, SSC-W

</td>

</tr>

</tbody>

</table>

## Remove `H/W` channels

``` r
res <- cqc_match(groups, ref = 4) 
res
```

``` r
cqc_fix(res)
```

## Coerce it directly into single `GatingSet` (zero-copying)

``` r
gs <- merge_list_to_gs(cqc_data)
gs
```

    ## A GatingSet with 63 samples
