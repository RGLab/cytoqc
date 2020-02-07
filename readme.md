
# cytoqc – A standardization tool for openCyto

*cytoqc* checks and standardizes channels, markers, keywords, gates of
the cytodata .

``` r
library(flowCore)
library(flowWorkspace)
library(cytoqc)
```

## Load the FCS

``` r
files <- list.files(data_dir, ".fcs", full.names = TRUE)
cqc_data <- cqc_load_fcs(files)
cqc_data
```

    ## cytoqc data: 
    ## 21  samples

## The basic workflow can be summarised as four steps:

1.  check
2.  match
3.  fix

### 1\. Check the consistency across samples

``` r
check_results <- cqc_check(cqc_data, type = "channel")
check_results
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nFCS

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

channel

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:left;">

FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, FSC-H, SSC-H, Time

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

channelA, FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, FSC-H, SSC1-H, Time

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

FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, fsc-h, SSC-H

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

FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, fsc-h, SSC1-H, Time

</td>

</tr>

</tbody>

</table>

### 2\. Match the reference

``` r
res <- cqc_match(check_results, ref = 3) 
res
```

    ##                Ref        1     2      4
    ## 1            FL1-H        ✓     ✓      ✓
    ## 2            FL2-A        ✓     ✓      ✓
    ## 3            FL2-H        ✓     ✓      ✓
    ## 4            FL3-H        ✓     ✓      ✓
    ## 5            FL4-H        ✓     ✓      ✓
    ## 6            FSC-H        ✓ fsc-h  fsc-h
    ## 7            SSC-H   SSC1-H     ✓ SSC1-H
    ## 8             Time        ✓  <NA>      ✓
    ## 10 To Delete       channelA

### 3\. Apply the fix

``` r
cqc_fix(res)
```

## Update check report

``` r
check_results <- cqc_check(cqc_data, type = "channel")
check_results
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nFCS

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

channel

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:left;">

FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, FSC-H, SSC-H, Time

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

FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, FSC-H, SSC-H

</td>

</tr>

</tbody>

</table>

## Drop outlier group

``` r
check_results <- cqc_drop_groups(check_results, id = 1)
check_results
```

<table class="table table-bordered" style="font-size: 12px; width: auto !important; ">

<thead>

<tr>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

group\_id

</th>

<th style="text-align:right;color: black !important;background-color: #e5f5e0 !important;">

nFCS

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

channel

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:left;">

FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, FSC-H, SSC-H, Time

</td>

</tr>

</tbody>

</table>

## Return the cleaned data

``` r
cqc_data <- cqc_get_data(check_results)
cqc_data
```

    ## cytoqc data: 
    ## 20  samples

## Coerce it inot `cytoset`

``` r
cytoset(cqc_data)
```

    ## A cytoset with 20 samples.
    ## 
    ##   column names:
    ##     FSC-H, SSC-H, FL1-H, FL2-H, FL3-H, FL2-A, FL4-H, Time
