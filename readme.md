
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
3.  recommend
4.  fix

### 1\. Check the consistency across samples

``` r
check_results <- cqc_check(cqc_data, type = "channel")
summary(check_results)
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

channelA,SSC1-H

</td>

<td style="text-align:left;">

SSC-H

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

2

</td>

<td style="text-align:left;">

fsc-h

</td>

<td style="text-align:left;">

FSC-H,Time

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

4

</td>

<td style="text-align:left;">

fsc-h,SSC1-H

</td>

<td style="text-align:left;">

FSC-H,SSC-H

</td>

</tr>

</tbody>

</table>

### 3\. Recommend the fix

``` r
solution <- cqc_recommend(res)
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

SSC1-H –\> SSC-H

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="     text-decoration: line-through;">channelA</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

fsc-h –\> FSC-H

</td>

</tr>

</tbody>

</table>

### 4\. Apply the fix

``` r
cqc_fix(solution)
```

## Update check report

``` r
check_results <- cqc_check(cqc_data, type = "channel")
summary(check_results)
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
summary(check_results)
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
