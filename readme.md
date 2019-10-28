
# cytoqc – A QC tool for openCyto

*cytoqc* performs pre-cleaning, pre-gating and post-gating QC checks.

``` r
library(flowCore)
library(flowWorkspace)
library(cytoqc)
```

## Load the FCS

``` r
files <- list.files(data_dir, ".fcs", full.names = TRUE)
cq_data <- cq_load_fcs(files)
cq_data
```

    ## cytoqc data: 
    ## 21  samples

## QC for channels

### Determin the reference channels

``` r
reference <- cq_find_reference_params(cq_data, type = "channel")
paste(reference, collapse = ", ")
```

    ## [1] "FL1-H, FL2-A, FL2-H, FL3-H, FL4-H, FSC-H, SSC-H, Time"

### QC check

``` r
check_results <- cq_check_params(cq_data, reference, type = "channel")
format(check_results)
```

<table class="table table-bordered" style="width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

FCS

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

unknown

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

missing

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

s6a01.fcs

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

s6a02.fcs

</td>

<td style="text-align:left;">

SSC1-H,channelA

</td>

<td style="text-align:left;">

SSC-H

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

s6a03.fcs

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

### Propose the fix

``` r
solution <- cq_fix_param_solution(check_results) 
format(solution)
```

<table class="table table-bordered" style="width: auto !important; ">

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

fsc-h –\> FSC-H

</td>

</tr>

<tr>

<td style="text-align:left;">

SSC1-H –\> SSC-H

</td>

</tr>

</tbody>

</table>

*Show the itemized
details*

``` r
format(solution, itemize = TRUE)
```

<table class="table table-bordered" style="width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

FCS

</th>

<th style="text-align:left;color: black !important;background-color: #e5f5e0 !important;">

Proposed change

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

s6a01.fcs

</td>

<td style="text-align:left;">

fsc-h –\> FSC-H

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

s6a02.fcs

</td>

<td style="text-align:left;">

SSC1-H –\>
SSC-H

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;font-weight: bold;" rowspan="2">

s6a03.fcs

</td>

<td style="text-align:left;">

fsc-h –\> FSC-H

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

SSC1-H –\> SSC-H

</td>

</tr>

</tbody>

</table>

### Apply the fix

After reviewing the `solution` (revise it if needed), pass it to the
`cq_fix_params`

``` r
cq_fix_params(cq_data, solution)
```

### Update QC report

``` r
check_results <- cq_check_params(cq_data, reference, type = "channel")
format(check_results)
```

<table class="table table-bordered" style="width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

FCS

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

unknown

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

missing

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

s6a01.fcs

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

Time

</td>

</tr>

<tr>

<td style="text-align:left;font-weight: bold;">

s6a02.fcs

</td>

<td style="text-align:left;">

channelA

</td>

<td style="text-align:left;">

</td>

</tr>

</tbody>

</table>

### Drop redundant channel

``` r
cq_data <- cq_drop_redundant_params(cq_data, check_results)
```

### Refresh QC report and drop the samples that still can not be fixed

``` r
check_results <- cq_check_params(cq_data, reference, type = "channel")
format(check_results)
```

<table class="table table-bordered" style="width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

FCS

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

unknown

</th>

<th style="text-align:left;color: black !important;background-color: #9ebcda !important;">

missing

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-weight: bold;">

s6a01.fcs

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

Time

</td>

</tr>

</tbody>

</table>

``` r
cq_data <- cq_drop_samples(cq_data, check_results)
length(cq_data)
```

    ## [1] 20

### QC for marker

``` r
reference <- cq_find_reference_params(cq_data, type = "marker")
check_results <- cq_check_params(cq_data, reference, type = "marker")
format(check_results)
```
