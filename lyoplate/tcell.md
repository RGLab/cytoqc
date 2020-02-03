
# QC for tcell panel

``` r
library(flowCore)
library(flowWorkspace)
library(cytoqc)
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
```

## QC Check gates

``` r
cqc_data <- cqc_gs_list(gslist)

#group by gates
groups <- cqc_check(cqc_data, "gate")
su <- summary(groups)
knit_print(su)
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
knit_print(diff(su))
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
plot_diff(groups)
```

![](tcell_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](tcell_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# match reference
match_result <- cqc_match(groups, ref = 1, select = c(2))
match_result
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

group\_id

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Not in
reference

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Missing channels

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

2

</td>

<td style="text-align:left;">

singlets

</td>

<td style="text-align:left;">

</td>

</tr>

</tbody>

</table>

``` r
solution <- cqc_recommend(match_result)
solution
```

<table class="table table-bordered table-condensed" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

Proposed change

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

<span style="     text-decoration: line-through;">singlets</span>

</td>

</tr>

</tbody>

</table>

``` r
cqc_fix(solution)

summary(cqc_check(cqc_data, "gate"))
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
su <- summary(groups)
su
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
\<V450-A\>, \<V545-A\>, FSC-A, FSC-H, FSC-W, SSC-A, SSC-H, SSC-W,
Time

</td>

</tr>

</tbody>

</table>

``` r
diff(su)
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
Blue-A\>, \<PE YG-A\>, \<PE-Cy7 YG-A\>, \<PerCP-Cy5-5-A\>

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
\<PE-A\>, \<PE-Cy7-A\>, \<PerCP-Cy5-5-A\>, FSC-H, FSC-W, SSC-H, SSC-W

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
\<PE Cy7 YG-A\>, \<PE-A\>, \<PerCP-Cy5-5-A\>, FSC-W, SSC-W

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
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>

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
\<FITC-A\>, \<PE-A\>, \<PE-Cy7-A\>, \<PerCP-Cy5-5-A\>, FSC-H, FSC-W,
SSC-H, SSC-W

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
\<PerCP-Cy5-5-A\>, \<V450-A\>, \<V500-A\>, FSC-H, SSC-H

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
\<V450-A\>, \<V545-A\>, FSC-H, FSC-W, SSC-H, SSC-W

</td>

</tr>

</tbody>

</table>

## Channels are very different across centers so move on to check marker

``` r
#summarise marker groups
groups <- cqc_check(cqc_data, "marker")
su <- summary(groups)
su
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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

CD197, CD3, CD38, CD4, CD45RA, CD8, HLA-DR, LIVE
DEAD

</td>

</tr>

</tbody>

</table>

``` r
knit_print(diff(su))
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
res <- cqc_match(groups, ref = 3, select = c(1,2,4,5,6)) 
res
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

group\_id

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Not in
reference

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Missing channels

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

1

</td>

<td style="text-align:left;">

CCR7 PE,CD3 V450,CD38 APC,CD4 PerCP-Cy55,CD45RA PE-Cy7,CD8 APC-H7,HLA-DR
V500,Live Dead FITC

</td>

<td style="text-align:left;">

CCR7,CD3,CD38,CD4,CD45RA,CD8,HLA-DR,LIVE

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

2

</td>

<td style="text-align:left;">

HLA DR,Live/Dead

</td>

<td style="text-align:left;">

HLA-DR,LIVE

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

4

</td>

<td style="text-align:left;">

HLADR,LIVE\_GREEN

</td>

<td style="text-align:left;">

HLA-DR,LIVE

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

5

</td>

<td style="text-align:left;">

CD197,LIVE DEAD

</td>

<td style="text-align:left;">

CCR7,LIVE

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

6

</td>

<td style="text-align:left;">

CD197,Live Green

</td>

<td style="text-align:left;">

CCR7,LIVE

</td>

</tr>

</tbody>

</table>

``` r
solution <- cqc_recommend(res, max.distance = 0.6)
solution
```

<table class="table table-bordered table-condensed" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

Proposed change

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

CCR7 PE –\> CCR7

</td>

</tr>

<tr>

<td style="text-align:left;">

CD38 APC –\> CD38

</td>

</tr>

<tr>

<td style="text-align:left;">

CD3 V450 –\> CD3

</td>

</tr>

<tr>

<td style="text-align:left;">

HLA-DR V500 –\> HLA-DR

</td>

</tr>

<tr>

<td style="text-align:left;">

CD45RA PE-Cy7 –\> CD45RA

</td>

</tr>

<tr>

<td style="text-align:left;">

HLA DR –\> HLA-DR

</td>

</tr>

<tr>

<td style="text-align:left;">

Live/Dead –\> LIVE

</td>

</tr>

<tr>

<td style="text-align:left;">

HLADR –\> HLA-DR

</td>

</tr>

<tr>

<td style="text-align:left;">

LIVE\_GREEN –\> LIVE

</td>

</tr>

<tr>

<td style="text-align:left;">

CD197 –\> CCR7

</td>

</tr>

<tr>

<td style="text-align:left;">

LIVE DEAD –\> LIVE

</td>

</tr>

<tr>

<td style="text-align:left;">

Live Green –\> LIVE

</td>

</tr>

</tbody>

</table>

``` r
cqc_fix(solution)
```

## update checks

``` r
groups <- cqc_check(cqc_data, "marker")
summary(groups)
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

marker

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CCR7, CD3, CD38, CD4, CD45RA, CD8, HLA-DR, LIVE

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

CCR7, CD3, CD38, CD4 PerCP-Cy55, CD45RA, CD8 APC-H7, HLA-DR, Live Dead
FITC

</td>

</tr>

</tbody>

</table>

## Relax string matching to clean up the rest

``` r
res <- cqc_match(groups, ref = 2, select = c(1)) 
res
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

group\_id

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Not in
reference

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Missing channels

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

1

</td>

<td style="text-align:left;">

CD4 PerCP-Cy55,CD8 APC-H7,Live Dead FITC

</td>

<td style="text-align:left;">

CD4,CD8,LIVE

</td>

</tr>

</tbody>

</table>

``` r
solution <- cqc_recommend(res, max.distance = 0.8)
solution
```

<table class="table table-bordered table-condensed" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

Proposed change

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

CD8 APC-H7 –\> CD8

</td>

</tr>

<tr>

<td style="text-align:left;">

Live Dead FITC –\> LIVE

</td>

</tr>

<tr>

<td style="text-align:left;">

CD4 PerCP-Cy55 –\> CD4

</td>

</tr>

</tbody>

</table>

``` r
cqc_fix(solution)

groups <- cqc_check(cqc_data, "marker")
summary(groups)
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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

CCR7, CD3, CD38, CD4, CD45RA, CD8, HLA-DR,
LIVE

</td>

</tr>

</tbody>

</table>

## Use the marker as reference to standardize the channels across centers

``` r
panel <- cf_get_panel(gs_cyto_data(cqc_data[[1]])[[1]], skip_na = TRUE)
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
summary(groups)
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

channel

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

marker

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

\<APC-A\>

</td>

<td style="text-align:left;">

CD38

</td>

<td style="text-align:right;vertical-align: top !important;" rowspan="8">

1

</td>

<td style="text-align:right;vertical-align: top !important;" rowspan="8">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

\<APC-H7-A\>

</td>

<td style="text-align:left;">

CD8

</td>

</tr>

<tr>

<td style="text-align:left;">

\<FITC-A\>

</td>

<td style="text-align:left;">

LIVE

</td>

</tr>

<tr>

<td style="text-align:left;">

\<PE-A\>

</td>

<td style="text-align:left;">

CCR7

</td>

</tr>

<tr>

<td style="text-align:left;">

\<PE-Cy7-A\>

</td>

<td style="text-align:left;">

CD45RA

</td>

</tr>

<tr>

<td style="text-align:left;">

\<PerCP-Cy5-5-A\>

</td>

<td style="text-align:left;">

CD4

</td>

</tr>

<tr>

<td style="text-align:left;">

\<V450-A\>

</td>

<td style="text-align:left;">

CD3

</td>

</tr>

<tr>

<td style="text-align:left;">

\<V500-A\>

</td>

<td style="text-align:left;">

HLA-DR

</td>

</tr>

</tbody>

</table>

## Refresh QC report

``` r
groups <- cqc_check(cqc_data, "channel")
su <- summary(groups)
su
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
knit_print(diff(su))
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nObject

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
res <- cqc_match(groups, ref = 4, select = c(1,2,3)) 
res
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

group\_id

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Not in
reference

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

Missing channels

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

1

</td>

<td style="text-align:left;">

FSC-H,FSC-W,SSC-H,SSC-W

</td>

<td style="text-align:left;">

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

2

</td>

<td style="text-align:left;">

FSC-H,SSC-H

</td>

<td style="text-align:left;">

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

3

</td>

<td style="text-align:left;">

FSC-W,SSC-W

</td>

<td style="text-align:left;">

</td>

</tr>

</tbody>

</table>

``` r
solution <- cqc_recommend(res, max.distance = 0.6)
solution
```

<table class="table table-bordered table-condensed" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

Proposed change

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

<span style="     text-decoration: line-through;">FSC-H</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="     text-decoration: line-through;">FSC-W</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="     text-decoration: line-through;">SSC-H</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="     text-decoration: line-through;">SSC-W</span>

</td>

</tr>

</tbody>

</table>

``` r
cqc_fix(solution)
```

## Coerce it directly into single `GatingSet` (zero-copying)

``` r
gs <- merge_list_to_gs(cqc_data)
gs
```

    ## A GatingSet with 63 samples
