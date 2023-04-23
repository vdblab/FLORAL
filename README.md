
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FLORAL: Fit LOg-RAtio Lasso regression for compositional covariates

<!-- badges: start -->

[![R-CMD-check](https://github.com/vdblab/FLORAL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vdblab/FLORAL/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/vdblab/FLORAL/branch/master/graph/badge.svg)](https://app.codecov.io/gh/vdblab/FLORAL?branch=master)
<!-- badges: end -->

`FLORAL` is an open-source computational tool to perform log-ratio lasso
regression modeling and microbial feature selection for continuous,
binary, time-to-event, and competing risk outcomes. The proposed method
adapts the augmented Lagrangian algorithm for a zero-sum constraint
optimization problem while enabling a two-stage screening process for
extended false-positive control.

## Installation

You can install the development version of `FLORAL` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vdblab/FLORAL")
```

## Example

Here is a toy example for linear regression with 10-fold
cross-validation for a simulated data with 50 samples and 100
compositional features. Option `progress=TRUE` can be used to show the
progress bar of the running algorithm.

The data simulation procedure is described in the preprint.

``` r
set.seed(23420)
library(FLORAL)

dat <- simu(n=50,p=100,model="linear")
fit <- FLORAL(dat$xcount,dat$y,family="gaussian",ncv=10,progress=FALSE,table=TRUE)
```

To view plots of cross-validated prediction error and parameter
coefficients, use `fit$pmse` or `fit$pcoef`:

<img src="man/figures/README-plot-1.png" width="100%" /><img src="man/figures/README-plot-2.png" width="100%" />

To view selected compositional features, use `fit$selected.feature`,
where features are sorted by their names. Features under `min` and `1se`
correspond to penalty parameter *λ*<sub>min</sub> and *λ*<sub>1se</sub>,
respectively. Features under `min.2stage` and `1se.2stage` are obtained
after applying 2-stage filtering based on features under `min` and
`1se`, respectively.

We recommend interpreting the selected compositional features as
potential predictive markers to the outcome in the regression model in
the sense that the cross-validated prediction error is improved by
considering these selected features.

``` r
fit$selected.feature
#> $min
#>  [1] "taxa1"  "taxa10" "taxa13" "taxa2"  "taxa20" "taxa3"  "taxa32" "taxa39"
#>  [9] "taxa5"  "taxa6"  "taxa60" "taxa7"  "taxa75" "taxa76" "taxa79" "taxa8" 
#> [17] "taxa84" "taxa9"  "taxa92"
#> 
#> $`1se`
#>  [1] "taxa1"  "taxa10" "taxa13" "taxa2"  "taxa20" "taxa3"  "taxa32" "taxa39"
#>  [9] "taxa5"  "taxa6"  "taxa7"  "taxa75" "taxa8"  "taxa84" "taxa9" 
#> 
#> $min.2stage
#>  [1] "taxa1"  "taxa10" "taxa13" "taxa2"  "taxa20" "taxa3"  "taxa32" "taxa5" 
#>  [9] "taxa6"  "taxa60" "taxa7"  "taxa79" "taxa8"  "taxa84" "taxa9"  "taxa92"
#> 
#> $`1se.2stage`
#>  [1] "taxa1"  "taxa10" "taxa13" "taxa2"  "taxa20" "taxa3"  "taxa32" "taxa5" 
#>  [9] "taxa6"  "taxa7"  "taxa8"  "taxa84" "taxa9"
```

To get specific log-ratios selected by the 2-stage procedure, use
`fit$step2.log-ratios`, where `min` and `1se` display the log-ratios
between features. For each identified ratio, `min.idx` and `1se.idx`
return the column indices in the original input matrix for the two
corresponding features forming the ratio.

``` r
fit$step2.ratios
#> $min
#>  [1] "taxa1/taxa13" "taxa1/taxa20" "taxa1/taxa84" "taxa2/taxa5"  "taxa3/taxa8" 
#>  [6] "taxa3/taxa92" "taxa5/taxa8"  "taxa6/taxa9"  "taxa7/taxa10" "taxa7/taxa79"
#> [11] "taxa8/taxa60" "taxa9/taxa32" "taxa9/taxa92"
#> 
#> $`1se`
#>  [1] "taxa1/taxa13" "taxa1/taxa20" "taxa1/taxa84" "taxa2/taxa5"  "taxa3/taxa8" 
#>  [6] "taxa5/taxa8"  "taxa6/taxa7"  "taxa6/taxa9"  "taxa7/taxa10" "taxa9/taxa32"
#> 
#> $min.idx
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> [1,]    1    1    1    2    3    3    5    6    7     7     8     9     9
#> [2,]   13   20   84    5    8   92    8    9   10    79    60    32    92
#> 
#> $`1se.idx`
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    1    1    2    3    5    6    6    7     9
#> [2,]   13   20   84    5    8    8    7    9   10    32
```

More detailed interpretations can be obtained for the selected
log-ratios. First, the selected log-ratios also improve the
cross-validated prediction errors because these log-ratios are derived
from the constrained lasso estimate. Moreover, as guided by the
association table between log-ratios and the outcome, it is possible to
interpret the directions of the covariate effects associated with
certain log-ratios on the outcome. To view detailed associations between
selected log-ratios and the outcome, set `table=TRUE` and use
`fit$step2.tables` to print `gtsummary` tables for the multivariable
stepwise regression models obtained by the 2-stage procedure.

``` r
fit$step2.tables$min
```

<div id="tqlkaqhafh" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#tqlkaqhafh .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#tqlkaqhafh .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#tqlkaqhafh .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#tqlkaqhafh .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#tqlkaqhafh .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tqlkaqhafh .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#tqlkaqhafh .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#tqlkaqhafh .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#tqlkaqhafh .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#tqlkaqhafh .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#tqlkaqhafh .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#tqlkaqhafh .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#tqlkaqhafh .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#tqlkaqhafh .gt_from_md > :first-child {
  margin-top: 0;
}

#tqlkaqhafh .gt_from_md > :last-child {
  margin-bottom: 0;
}

#tqlkaqhafh .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#tqlkaqhafh .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#tqlkaqhafh .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#tqlkaqhafh .gt_row_group_first td {
  border-top-width: 2px;
}

#tqlkaqhafh .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tqlkaqhafh .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#tqlkaqhafh .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#tqlkaqhafh .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tqlkaqhafh .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tqlkaqhafh .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#tqlkaqhafh .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#tqlkaqhafh .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tqlkaqhafh .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#tqlkaqhafh .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#tqlkaqhafh .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#tqlkaqhafh .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#tqlkaqhafh .gt_left {
  text-align: left;
}

#tqlkaqhafh .gt_center {
  text-align: center;
}

#tqlkaqhafh .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#tqlkaqhafh .gt_font_normal {
  font-weight: normal;
}

#tqlkaqhafh .gt_font_bold {
  font-weight: bold;
}

#tqlkaqhafh .gt_font_italic {
  font-style: italic;
}

#tqlkaqhafh .gt_super {
  font-size: 65%;
}

#tqlkaqhafh .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#tqlkaqhafh .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#tqlkaqhafh .gt_indent_1 {
  text-indent: 5px;
}

#tqlkaqhafh .gt_indent_2 {
  text-indent: 10px;
}

#tqlkaqhafh .gt_indent_3 {
  text-indent: 15px;
}

#tqlkaqhafh .gt_indent_4 {
  text-indent: 20px;
}

#tqlkaqhafh .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col"><strong>Characteristic</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>Beta</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>95% CI</strong><sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>p-value</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">taxa1/taxa13</td>
<td class="gt_row gt_center">0.09</td>
<td class="gt_row gt_center">0.05, 0.14</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa1/taxa20</td>
<td class="gt_row gt_center">0.06</td>
<td class="gt_row gt_center">0.00, 0.11</td>
<td class="gt_row gt_center">0.040</td></tr>
    <tr><td class="gt_row gt_left">taxa1/taxa84</td>
<td class="gt_row gt_center">0.05</td>
<td class="gt_row gt_center">0.00, 0.10</td>
<td class="gt_row gt_center">0.064</td></tr>
    <tr><td class="gt_row gt_left">taxa2/taxa5</td>
<td class="gt_row gt_center">-0.10</td>
<td class="gt_row gt_center">-0.13, -0.06</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa3/taxa8</td>
<td class="gt_row gt_center">0.06</td>
<td class="gt_row gt_center">0.01, 0.11</td>
<td class="gt_row gt_center">0.019</td></tr>
    <tr><td class="gt_row gt_left">taxa3/taxa92</td>
<td class="gt_row gt_center">0.07</td>
<td class="gt_row gt_center">0.01, 0.12</td>
<td class="gt_row gt_center">0.021</td></tr>
    <tr><td class="gt_row gt_left">taxa5/taxa8</td>
<td class="gt_row gt_center">0.10</td>
<td class="gt_row gt_center">0.06, 0.15</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa6/taxa9</td>
<td class="gt_row gt_center">-0.16</td>
<td class="gt_row gt_center">-0.20, -0.13</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa7/taxa10</td>
<td class="gt_row gt_center">0.17</td>
<td class="gt_row gt_center">0.12, 0.22</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa7/taxa79</td>
<td class="gt_row gt_center">0.05</td>
<td class="gt_row gt_center">0.00, 0.10</td>
<td class="gt_row gt_center">0.083</td></tr>
    <tr><td class="gt_row gt_left">taxa8/taxa60</td>
<td class="gt_row gt_center">-0.06</td>
<td class="gt_row gt_center">-0.11, -0.01</td>
<td class="gt_row gt_center">0.021</td></tr>
    <tr><td class="gt_row gt_left">taxa9/taxa32</td>
<td class="gt_row gt_center">0.08</td>
<td class="gt_row gt_center">0.03, 0.14</td>
<td class="gt_row gt_center">0.005</td></tr>
    <tr><td class="gt_row gt_left">taxa9/taxa92</td>
<td class="gt_row gt_center">-0.03</td>
<td class="gt_row gt_center">-0.08, 0.02</td>
<td class="gt_row gt_center">0.2</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="4"><sup class="gt_footnote_marks">1</sup> CI = Confidence Interval</td>
    </tr>
  </tfoot>
</table>
</div>

``` r
fit$step2.tables$`1se`
```

<div id="mjbnlivwgz" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#mjbnlivwgz .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#mjbnlivwgz .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mjbnlivwgz .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#mjbnlivwgz .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#mjbnlivwgz .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mjbnlivwgz .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mjbnlivwgz .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#mjbnlivwgz .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#mjbnlivwgz .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#mjbnlivwgz .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#mjbnlivwgz .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#mjbnlivwgz .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#mjbnlivwgz .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#mjbnlivwgz .gt_from_md > :first-child {
  margin-top: 0;
}

#mjbnlivwgz .gt_from_md > :last-child {
  margin-bottom: 0;
}

#mjbnlivwgz .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#mjbnlivwgz .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#mjbnlivwgz .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#mjbnlivwgz .gt_row_group_first td {
  border-top-width: 2px;
}

#mjbnlivwgz .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mjbnlivwgz .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#mjbnlivwgz .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#mjbnlivwgz .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mjbnlivwgz .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mjbnlivwgz .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#mjbnlivwgz .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#mjbnlivwgz .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mjbnlivwgz .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mjbnlivwgz .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mjbnlivwgz .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mjbnlivwgz .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mjbnlivwgz .gt_left {
  text-align: left;
}

#mjbnlivwgz .gt_center {
  text-align: center;
}

#mjbnlivwgz .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#mjbnlivwgz .gt_font_normal {
  font-weight: normal;
}

#mjbnlivwgz .gt_font_bold {
  font-weight: bold;
}

#mjbnlivwgz .gt_font_italic {
  font-style: italic;
}

#mjbnlivwgz .gt_super {
  font-size: 65%;
}

#mjbnlivwgz .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#mjbnlivwgz .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#mjbnlivwgz .gt_indent_1 {
  text-indent: 5px;
}

#mjbnlivwgz .gt_indent_2 {
  text-indent: 10px;
}

#mjbnlivwgz .gt_indent_3 {
  text-indent: 15px;
}

#mjbnlivwgz .gt_indent_4 {
  text-indent: 20px;
}

#mjbnlivwgz .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col"><strong>Characteristic</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>Beta</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>95% CI</strong><sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>p-value</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">taxa1/taxa13</td>
<td class="gt_row gt_center">0.07</td>
<td class="gt_row gt_center">0.02, 0.12</td>
<td class="gt_row gt_center">0.007</td></tr>
    <tr><td class="gt_row gt_left">taxa1/taxa20</td>
<td class="gt_row gt_center">0.08</td>
<td class="gt_row gt_center">0.02, 0.13</td>
<td class="gt_row gt_center">0.007</td></tr>
    <tr><td class="gt_row gt_left">taxa1/taxa84</td>
<td class="gt_row gt_center">0.05</td>
<td class="gt_row gt_center">0.00, 0.10</td>
<td class="gt_row gt_center">0.069</td></tr>
    <tr><td class="gt_row gt_left">taxa2/taxa5</td>
<td class="gt_row gt_center">-0.09</td>
<td class="gt_row gt_center">-0.12, -0.06</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa3/taxa8</td>
<td class="gt_row gt_center">0.11</td>
<td class="gt_row gt_center">0.08, 0.14</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa5/taxa8</td>
<td class="gt_row gt_center">0.11</td>
<td class="gt_row gt_center">0.07, 0.15</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa6/taxa7</td>
<td class="gt_row gt_center">-0.05</td>
<td class="gt_row gt_center">-0.11, 0.00</td>
<td class="gt_row gt_center">0.066</td></tr>
    <tr><td class="gt_row gt_left">taxa6/taxa9</td>
<td class="gt_row gt_center">-0.13</td>
<td class="gt_row gt_center">-0.18, -0.07</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa7/taxa10</td>
<td class="gt_row gt_center">0.17</td>
<td class="gt_row gt_center">0.12, 0.22</td>
<td class="gt_row gt_center"><0.001</td></tr>
    <tr><td class="gt_row gt_left">taxa9/taxa32</td>
<td class="gt_row gt_center">0.09</td>
<td class="gt_row gt_center">0.03, 0.15</td>
<td class="gt_row gt_center">0.006</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="4"><sup class="gt_footnote_marks">1</sup> CI = Confidence Interval</td>
    </tr>
  </tfoot>
</table>
</div>

For binary and survival outcomes, please specify `family="binomial"` and
`family="cox"` accordingly.

``` r
dat.bin <- simu(n=50,p=100,model="binomial")
fit.bin <- FLORAL(dat.bin$xcount,dat.bin$y,family="binomial",ncv=10,progress=FALSE)

dat.cox <- simu(n=50,p=100,model="cox")
fit.cox <- FLORAL(dat.cox$xcount,survival::Surv(dat.cox$t,dat.cox$d),family="cox",ncv=10,progress=FALSE)
```
