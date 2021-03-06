---
title: "RMDdraft1MLMS"
output:
  minidown::mini_document:
    framework: water
    theme: light
    code_folding:
      source: show
      output: show
      message: hide
      warning: hide
      error: show
    toc: true
    toc_float: true
    toc_highlight: true
    tabset: true
    number_sections: true
    anchor_sections: false
    self_contained: false
    code_download: true
    math: "katex"
    keep_md: true
    fig_caption: true
vignette: >
  %\VignetteIndexEntry{Writing Vignette with the 'minidown' Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
---




This article introduces the MLMS R library and is a draft of a vignette[@R-base].

The github link for MLMS is 
[github.com/lilyacb/MLMS](https://github.com/lilyacb/MLMS).  Open this document in a browser for the links to work.

# Accessing information and making tables

Tables summarizing .dxf file contents can be produced with *read_summary()*.

```{.r .details .show summary='Source'}
  datafile<-"170525_NaHCO3 L + NaCl L_.dxf"
  #data_file<-"https://github.com/lilyacb/MLMS/blob/main/Data/170525_NaHCO3%2#0L%20%2B%20NaCl%20L_.dxf"  # would work if changed to csv!!
  #file.summ<-read_csv(data_file)
  # Print table without kable
  fileI.summ<-read_summary(datafile)
```

```{.details .show summary='Output'}
##                   Length Class           Mode
## version            1     package_version list
## read_options       4     -none-          list
## file_info         16     tbl_df          list
## method_info        3     -none-          list
## raw_data           5     tbl_df          list
## vendor_data_table 39     tbl_df          list
```
<br><br>

You can use kable to make a fancier table.

```{.r .details .show summary='Source'}
  #class.source='details hide'
 # Using kable
  knitr::kable(fileI.summ,caption="170525_NaHCO3 L + NaCl L_.dxf file summary")
```

<table>
<caption>170525_NaHCO3 L + NaCl L_.dxf file summary</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Length </th>
   <th style="text-align:left;"> Class </th>
   <th style="text-align:left;"> Mode </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> version </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> package_version </td>
   <td style="text-align:left;"> list </td>
  </tr>
  <tr>
   <td style="text-align:left;"> read_options </td>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> -none- </td>
   <td style="text-align:left;"> list </td>
  </tr>
  <tr>
   <td style="text-align:left;"> file_info </td>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:left;"> tbl_df </td>
   <td style="text-align:left;"> list </td>
  </tr>
  <tr>
   <td style="text-align:left;"> method_info </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> -none- </td>
   <td style="text-align:left;"> list </td>
  </tr>
  <tr>
   <td style="text-align:left;"> raw_data </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> tbl_df </td>
   <td style="text-align:left;"> list </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vendor_data_table </td>
   <td style="text-align:left;"> 39 </td>
   <td style="text-align:left;"> tbl_df </td>
   <td style="text-align:left;"> list </td>
  </tr>
</tbody>
</table>
<br><br>

You can print information contained in the file_info, vendor_info and raw_data tabs of the .dxf file.

Get file information with *file_info()*.

```{.r .details .show summary='Source'}
# Can get file information
fi.df<-file_info(files=datafile)
knitr::kable(fi.df,caption="170525_NaHCO3 L + NaCl L_.dxf file information")
```

<table>
<caption>170525_NaHCO3 L + NaCl L_.dxf file information</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> file_id </th>
   <th style="text-align:left;"> Identifier_1 </th>
   <th style="text-align:left;"> Analysis </th>
   <th style="text-align:left;"> Preparation </th>
   <th style="text-align:left;"> Date_and_Time </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:left;"> NaHCO3 L + NaCl L </td>
   <td style="text-align:left;"> 4172 </td>
   <td style="text-align:left;"> 2% CO2 in He 24hrs </td>
   <td style="text-align:left;"> 2017-05-25 21:15:20 </td>
  </tr>
</tbody>
</table>
<br>

Get the vendor data table with *vendor_info()*.

```{.r .details .show summary='Source'}
# Can get vendor info
vi.df<-vendor_info(datafile)
kbl(head(vi.df)[1:3,],caption="170525_NaHCO3 L + NaCl L_.dxf vendor data") %>% 
  kable_paper() %>%
  scroll_box(width=5,height = "200px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:200px; overflow-x: scroll; width:5; "><table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
<caption>170525_NaHCO3 L + NaCl L_.dxf vendor data</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Identifier_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Peak_Nr </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Start </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Rt </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> End </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Intensity_All </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> rIntensity_All </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Ampl_44 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Ampl_45 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Ampl_46 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> d13C/12C </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> d18O/16O </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NaHCO3 L + NaCl L </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 27.1700000762939 </td>
   <td style="text-align:left;"> 47.443000793457 </td>
   <td style="text-align:left;"> 50.3689994812012 </td>
   <td style="text-align:left;"> 40.90863519699 </td>
   <td style="text-align:left;"> 141539.302508887 </td>
   <td style="text-align:left;"> 2062.13298836363 </td>
   <td style="text-align:left;"> 2371.49813464285 </td>
   <td style="text-align:left;"> 2812.49753618504 </td>
   <td style="text-align:left;"> -36.8309479493618 </td>
   <td style="text-align:left;"> -39.9994709523972 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NaHCO3 L + NaCl L </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> 67.088996887207 </td>
   <td style="text-align:left;"> 87.3619995117188 </td>
   <td style="text-align:left;"> 90.0790023803711 </td>
   <td style="text-align:left;"> 40.8018339597763 </td>
   <td style="text-align:left;"> 141166.762952385 </td>
   <td style="text-align:left;"> 2065.4354230274 </td>
   <td style="text-align:left;"> 2375.34645516822 </td>
   <td style="text-align:left;"> 2816.45019646327 </td>
   <td style="text-align:left;"> -36.9 </td>
   <td style="text-align:left;"> -40.0000000000001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NaHCO3 L + NaCl L </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 131.878997802734 </td>
   <td style="text-align:left;"> 134.177993774414 </td>
   <td style="text-align:left;"> 138.149002075195 </td>
   <td style="text-align:left;"> 0.775306018352309 </td>
   <td style="text-align:left;"> 2741.90298831581 </td>
   <td style="text-align:left;"> 497.765554776576 </td>
   <td style="text-align:left;"> 587.296395809005 </td>
   <td style="text-align:left;"> 684.745514792838 </td>
   <td style="text-align:left;"> -12.8581753517045 </td>
   <td style="text-align:left;"> -3.96450973206153 </td>
  </tr>
</tbody>
</table></div>
<br>

Get the raw data using *raw_data()*

```{.r .details .show summary='Source'}
# Can get the raw data
raw.df<-raw_data(datafile)
knitr::kable(head(raw.df),caption="170525_NaHCO3 L + NaCl L_.dxf raw data")
```

<table>
<caption>170525_NaHCO3 L + NaCl L_.dxf raw data</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> file_id </th>
   <th style="text-align:right;"> tp </th>
   <th style="text-align:right;"> time.s </th>
   <th style="text-align:right;"> v44.mV </th>
   <th style="text-align:right;"> v45.mV </th>
   <th style="text-align:right;"> v46.mV </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.209 </td>
   <td style="text-align:right;"> 1.412726 </td>
   <td style="text-align:right;"> 0.8119730 </td>
   <td style="text-align:right;"> 2.852197 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.418 </td>
   <td style="text-align:right;"> 1.416549 </td>
   <td style="text-align:right;"> 0.8024275 </td>
   <td style="text-align:right;"> 2.854114 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.627 </td>
   <td style="text-align:right;"> 1.414637 </td>
   <td style="text-align:right;"> 0.7928824 </td>
   <td style="text-align:right;"> 2.823437 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.836 </td>
   <td style="text-align:right;"> 1.412726 </td>
   <td style="text-align:right;"> 0.9628401 </td>
   <td style="text-align:right;"> 3.186850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1.045 </td>
   <td style="text-align:right;"> 1.406990 </td>
   <td style="text-align:right;"> 0.8768920 </td>
   <td style="text-align:right;"> 3.003719 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1.254 </td>
   <td style="text-align:right;"> 1.410814 </td>
   <td style="text-align:right;"> 0.7508881 </td>
   <td style="text-align:right;"> 2.712264 </td>
  </tr>
</tbody>
</table>
<br>

Get the resistor information using *resistor_data()*

```{.r .details .show summary='Source'}
# Can get the resistor information
resist<-resistor_data(datafile)
knitr::kable(resist,caption="170525_NaHCO3 L + NaCl L_.dxf resistor information")
```

<table>
<caption>170525_NaHCO3 L + NaCl L_.dxf resistor information</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> file_id </th>
   <th style="text-align:right;"> cup </th>
   <th style="text-align:right;"> R.Ohm </th>
   <th style="text-align:left;"> mass </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 3e+08 </td>
   <td style="text-align:left;"> 44 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3e+10 </td>
   <td style="text-align:left;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1e+11 </td>
   <td style="text-align:left;"> 46 </td>
  </tr>
</tbody>
</table>
<br>

Get the isotopic reference values with *reference_values_ratio()*

```{.r .details .show summary='Source'}
# Can get isotopic reference values with ratios
stand_ratio<-reference_values_ratio(datafile)
knitr::kable(stand_ratio,caption="170525_NaHCO3 L + NaCl L_.dxf isotopic reference values with ratios")
```

<table>
<caption>170525_NaHCO3 L + NaCl L_.dxf isotopic reference values with ratios</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> file_id </th>
   <th style="text-align:left;"> standard </th>
   <th style="text-align:left;"> gas </th>
   <th style="text-align:left;"> delta_name </th>
   <th style="text-align:right;"> delta_value </th>
   <th style="text-align:left;"> reference </th>
   <th style="text-align:left;"> element </th>
   <th style="text-align:left;"> ratio_name </th>
   <th style="text-align:right;"> ratio_value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:left;"> CO2_zero </td>
   <td style="text-align:left;"> CO2 </td>
   <td style="text-align:left;"> d 13C/12C </td>
   <td style="text-align:right;"> -36.9 </td>
   <td style="text-align:left;"> VPDB </td>
   <td style="text-align:left;"> C </td>
   <td style="text-align:left;"> R 13C/12C </td>
   <td style="text-align:right;"> 0.0111802 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:left;"> CO2_zero </td>
   <td style="text-align:left;"> CO2 </td>
   <td style="text-align:left;"> d 13C/12C </td>
   <td style="text-align:right;"> -36.9 </td>
   <td style="text-align:left;"> VPDB </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> R 18O/16O </td>
   <td style="text-align:right;"> 0.0020672 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:left;"> CO2_zero </td>
   <td style="text-align:left;"> CO2 </td>
   <td style="text-align:left;"> d 13C/12C </td>
   <td style="text-align:right;"> -36.9 </td>
   <td style="text-align:left;"> VPDB </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> R 17O/16O </td>
   <td style="text-align:right;"> 0.0003860 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:left;"> CO2_zero </td>
   <td style="text-align:left;"> CO2 </td>
   <td style="text-align:left;"> d 18O/16O </td>
   <td style="text-align:right;"> -40.0 </td>
   <td style="text-align:left;"> VSMOW </td>
   <td style="text-align:left;"> H </td>
   <td style="text-align:left;"> R 2H/1H </td>
   <td style="text-align:right;"> 0.0001558 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:left;"> CO2_zero </td>
   <td style="text-align:left;"> CO2 </td>
   <td style="text-align:left;"> d 18O/16O </td>
   <td style="text-align:right;"> -40.0 </td>
   <td style="text-align:left;"> VSMOW </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> R 17O/16O </td>
   <td style="text-align:right;"> 0.0003799 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 170525_NaHCO3 L + NaCl L_.dxf </td>
   <td style="text-align:left;"> CO2_zero </td>
   <td style="text-align:left;"> CO2 </td>
   <td style="text-align:left;"> d 18O/16O </td>
   <td style="text-align:right;"> -40.0 </td>
   <td style="text-align:left;"> VSMOW </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> R 18O/16O </td>
   <td style="text-align:right;"> 0.0020052 </td>
  </tr>
</tbody>
</table>
<br>

Use DT to render larger tables neatly. You can show only a few lines, have a search bar, filters and more.


```{.r .details .show summary='Source'}
datatable(vi.df,#filter="top",
          options=list(pageLength=5,scrollX=T))
```

```{=html}
<div id="htmlwidget-90e11e175380102ca2c9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-90e11e175380102ca2c9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"],["NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L","NaHCO3 L + NaCl L"],["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"],["27.1700000762939","67.088996887207","131.878997802734","166.572998046875","206.492004394531","271.282012939453","331.055999755859","390.621002197266","450.394989013672","510.169006347656","569.942993164062","629.716979980469","689.281982421875","749.056030273438","808.620971679688","843.315002441406"],["47.443000793457","87.3619995117188","134.177993774414","186.845993041992","226.555999755859","275.252990722656","334.817993164062","394.592010498047","454.365997314453","513.931030273438","573.705017089844","633.47900390625","693.252990722656","752.817993164062","812.591979980469","863.588012695312"],["50.3689994812012","90.0790023803711","138.149002075195","189.563003540039","229.481994628906","281.523010253906","341.088012695312","400.861999511719","460.427001953125","520.200988769531","579.765991210938","639.539978027344","699.10498046875","758.879028320312","818.443969726562","866.513977050781"],["40.90863519699","40.8018339597763","0.775306018352309","40.7667643372804","40.5869572281707","3.2905863791484","3.13872566031333","2.99676108613764","2.86171429394237","2.73679153407845","2.62090716786112","2.51637866908364","2.41483946341397","2.31754065522929","2.22933039823088","40.7124472200286"],["141539.302508887","141166.762952385","2741.90298831581","141046.972056521","140423.818600922","11634.7360679161","11097.4910169777","10595.5253245978","10116.641868754","9675.86233312315","9264.64028335832","8894.88525895595","8536.21830667717","8190.25073574662","7878.85016825028","140862.365302356"],["2062.13298836363","2065.4354230274","497.765554776576","2065.33132567091","2065.95164329982","835.629396537048","796.707468988884","760.583045951027","727.136230359232","695.731342346527","666.690982593185","639.816551703906","614.178994448332","590.535712752739","568.269862700224","2059.70138770054"],["2371.49813464285","2375.34645516822","587.296395809005","2375.28129849837","2375.8375421133","983.807398920264","938.38670479377","895.553766797341","855.914217894725","819.130433390809","784.737002483931","753.050001131472","722.521469497472","694.496962038834","668.386176906104","2368.64119184778"],["2812.49753618504","2816.45019646327","684.745514792838","2815.76340254216","2816.62504109653","1181.01548866047","1126.1267886412","1074.98412818962","1027.52041625678","983.555646869736","941.803287231497","904.883046181393","868.365553263493","833.933201559691","803.193882381927","2808.72003462978"],["-36.8309479493618","-36.9","-12.8581753517045","-36.8501122637148","-36.8704263734039","-13.3314973496644","-13.4239050632513","-13.6527361443061","-13.7448311807784","-13.798912139493","-13.889804179784","-14.0983142851769","-14.2407335597668","-14.5823805444807","-14.4921650499025","-36.817546278087"],["-39.9994709523972","-40.0000000000001","-3.96450973206153","-40.0105813222899","-40.0135757068457","-4.15988057764383","-4.15856399628201","-3.99410487962526","-4.27067517196855","-4.01350290214231","-4.3526332497722","-4.26916577111913","-4.08914826967777","-4.44953861274611","-4.4104799888387","-39.9775032876957"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Identifier_1<\/th>\n      <th>Peak_Nr<\/th>\n      <th>Start<\/th>\n      <th>Rt<\/th>\n      <th>End<\/th>\n      <th>Intensity_All<\/th>\n      <th>rIntensity_All<\/th>\n      <th>Ampl_44<\/th>\n      <th>Ampl_45<\/th>\n      <th>Ampl_46<\/th>\n      <th>d13C/12C<\/th>\n      <th>d18O/16O<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":5,"scrollX":true,"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}],"lengthMenu":[5,10,25,50,100]}},"evals":[],"jsHooks":[]}</script>
```
<br><br>

# MS data processing

## Peak areas

Peak areas were calculated using the *trapz* package, which implements numerical integration via the trapezoid rule.

Math with KaTex produces nice equations. You can highlight formulas or text.<br>
<br>
The trapezoid rule:
<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 20px;}
</style>
<div class = "blue">
$$
\int\limits_a^b f(x)dx\approx
\sum_{n=0}^{N-1}{\frac{1}{2}(f_n+f_{n+1})(\varDelta x)_n}
$$
</div>
<br>
You can insert images into an RMD doc.

```{.r .details .show summary='Source'}
trapImg<-readPNG("trapezoidalRuleImg.png")
grid.raster(trapImg)
```

![Graphical illustration of the trapezoid rule for numerical integration](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-9-1.png)
<br>

Calculate the peak areas using *trap_area_allPks()*.

```{.r .details .show summary='Source'}
# Can get peak areas
rawN<-"v44.mV"
areaPks<-trap_area_allPks(raw.df,vi.df,rawN)
knitr::kable(areaPks,caption="170525_NaHCO3 L + NaCl L_.dxf all peak areas")
```

<table>
<caption>170525_NaHCO3 L + NaCl L_.dxf all peak areas</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> Pk_Nr </th>
   <th style="text-align:right;"> trap_area </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 40.6839545 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 40.5842293 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.7984758 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 40.5428463 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 40.3672183 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3.3203306 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 3.1683712 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 3.0249778 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 2.8898409 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 2.7635269 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 2.6484939 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2.5426593 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 2.4407975 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 2.3417151 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 2.2528789 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 40.4814526 </td>
  </tr>
</tbody>
</table>
<br><br>


## Data visualization
<style>
p.caption {
  font-size: 0.9em;
  font-style: italic;
  color: grey;
  margin-right: 10%;
  margin-left: 10%;  
  text-align: justify;
}
</style>

You can graph the intensity vs retention time using *plot_ms()*.

```{.r .details .show summary='Source'}
# Intensity_All vs Rt for an experiment
ia.xname<-"Rt"
ia.yname<-"Intensity_All"
# plot
plot_ms(vi.df,ia.xname,ia.yname)
```

![Intensity_all vs Rt for 170525_NaHCO3_L_+_NaCl](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-11-1.png)
<br><br>


Plot individual peaks in an experiment with *plot_individual_peaks()*.

```{.r .details .show summary='Source'}
# Can plot peaks in an experiment individually (Intensity (mV) vs Rt)
time.s<-as.numeric(raw.df$time.s)
start.v1<-as.numeric(vi.df$Start)
end.v1<-as.numeric(vi.df$End)
v44<-as.numeric(raw.df$v44.mV)
# plot just the first peak to inspect
peak1.p<-plot_individual_peaks(start.v1,end.v1,time.s,v44,"1","v44.mV")
```

![Peak 1 intensity (v44.mV) vs. time](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-12-1.png)
<br><br>

Plot the raw data, intensity vs retention time, with *gg_raw_plot()*--make function!!

Use ggplot to plot the raw (**redo colour label in legend**)

```{.r .details .show summary='Source'}
raw.dat<-read.table("LLrawdat")
ggplot(raw.dat,aes(x=time.s,y=v44.mV))+
  geom_line(aes(color="v44.mV"))+
  geom_line(aes(x=time.s,y=v45.mV,color="v45.mV"))+
  geom_line(aes(x=time.s,y=v46.mV,color="v46.mV"))+
  labs(title="170525_NaHCO3 L + NaCl U",x="time.s",y="v44_v45_v46.mV")
```

![](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-13-1.png)<!-- -->
<br>

Plot the raw data of all files in a directory and export to a pdf file using *generic_plot_all_raw()* **change to ggplot version**--new func!

```{.r .details .show summary='Source'}
# Can plot raw data of all files in a directory
# Get all filenames of .dxf files in the directory
LLdir<-"NaHCO3_L_+_NaCl_L"
fileNames<-all_filenames(LLdir)
setwd(LLdir)
rawList<-raw_data_all(fileNames)
setwd("~/Desktop/EuropaMLMS/rmdMLMS")
# plot all raw data
generic_plot_all_raw(rawList)
```

![](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-14-1.png)<!-- -->![](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-14-2.png)<!-- -->![](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-14-3.png)<!-- -->![](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-14-4.png)<!-- -->![](MLMSprocessingTutorial_files/figure-html5/unnamed-chunk-14-5.png)<!-- -->
<br>

## Sorting and quality checks

### Automated sorting of a directory of .dxf files
A directory of .dxf files can be organized by Identifier_1 and Preparation method.
<br><br>



Contents of an unsorted directory of .dxf files.

```{.details .show summary='Output'}
## [1] "170525_NaHCO3 L + NaCl L_.dxf" "170525_NaHCO3 L + NaCl U_.dxf"
## [3] "170525_NaHCO3 L_.dxf"          "170525_NaHCO3 U + NaCl L_.dxf"
## [5] "170525_NaHCO3 U + NaCl U_.dxf" "170525_NaHCO3 U_.dxf"
```
<br>

Use *sort_by_identifier_1()* to sort a directory of .dxf files by Identifier_1

```{.r .details .show summary='Source'}
# Use Identifier_1 labels to sort data
unsortedPath<-"~/Desktop/EuropaMLMS/rmdMLMS/vignetteData/sortFolder"
setwd(unsortedPath)
sort_by_identifier_1(unsortedPath)
setwd("~/Desktop/EuropaMLMS/rmdMLMS")
```
<br>
Contents of the sorted directory.

```{.details .show summary='Output'}
## [1] "cache"             "NaHCO3 L"          "NaHCO3 L + NaCl L"
## [4] "NaHCO3 L + NaCl U" "NaHCO3 U"          "NaHCO3 U + NaCl L"
## [7] "NaHCO3 U + NaCl U"
```
<br>

### Automated quality-control and calibration of a directory of .dxf files

Perform quality checks and remove files that fail any of the checks.

<body>
**Quality checks:**
  <ul>
 1. 
 2. 
 3. 
 4. 
 5. 
 6.
 </ul>
<body>

The list above represents the order of quality checks and a file is removed immediately after it fails a check (before the next check is performed).

**1. Peaks present/number of peaks **<br>
Check that there are more than 0 and fewer than a specified number of peaks.


**2. Reference peaks**<br>
Check- 


**3. Reference peaks**<br>
checks-

<br>

#### Quality-control summary for a directory of .dxf files {.tabset .tabset-fade .tabset-pills}

Output from the quality check and calibration process can summarized in tables and analyzed. 

##### QC1 {.unnumbered}
(tab content)

##### QC2 {.unnumbered}
(tab content)

#### {.unnumbered}
Summaries for analyses that failed checks
<br><br>


#### Internal standards summary for a directory of .dxf files (if standards present)
Internal standards can be checked for a batch of runs (several .dxf files with internal standard experiments).
<br><br>


### Automated separation and quality-control of reference and sample peaks in a combined .csv file for multiple experiments
Processing of weekeqdata (three .csv files) into quality-checked reference and sample peaks 

<br><br>


# UMAP
Illustration of interactive graphics

## Unfiltered data
UMAP for different reactions that have not been quality-checked.<br>
This UMAP uses a directory of unsorted, unfiltered data and is for all the peaks in the data set.


```{.r .details .show summary='Source'}
unfiltPlot.dat<-read.table("UMAP_NaHCO3_+_NaCl_notQC")
p<-ggplot(unfiltPlot.dat,
       aes(x=x,y=y,color=Identifier_1))+
       geom_point()+
       labs(title="UMAP of NaHCO3 + NaCl (not quality-checked)")
ggplotly(p)
```

```{=html}
<div id="htmlwidget-32427431855900e79bb4" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-32427431855900e79bb4">{"x":{"data":[{"x":[10.3356594479537,10.2935309249608,-9.29826600296247,10.1284101282884,8.56211535341234,10.4261168175786,-9.39626719554616,-9.46177423465057,-9.30391095910075,-9.2031766597065,-9.48001412611369,-9.24662349387062,-9.44055292366261,-9.35950115742232,-9.27840942623549,7.59828571027462],"y":[3.8852469848589,3.95412847775669,0.197593842149332,3.82795658525943,3.98490394825349,4.17266369069777,0.790164448643495,0.643446148020958,0.303920437241228,0.260654392919486,0.036126558624446,-0.137981883983975,-0.26933954808535,-0.351811352951761,-0.303422016395703,5.44654491980404],"text":["x: 10.335659<br />y:   3.88524698<br />Identifier_1: NaHCO3 L","x: 10.293531<br />y:   3.95412848<br />Identifier_1: NaHCO3 L","x: -9.298266<br />y:   0.19759384<br />Identifier_1: NaHCO3 L","x: 10.128410<br />y:   3.82795659<br />Identifier_1: NaHCO3 L","x:  8.562115<br />y:   3.98490395<br />Identifier_1: NaHCO3 L","x: 10.426117<br />y:   4.17266369<br />Identifier_1: NaHCO3 L","x: -9.396267<br />y:   0.79016445<br />Identifier_1: NaHCO3 L","x: -9.461774<br />y:   0.64344615<br />Identifier_1: NaHCO3 L","x: -9.303911<br />y:   0.30392044<br />Identifier_1: NaHCO3 L","x: -9.203177<br />y:   0.26065439<br />Identifier_1: NaHCO3 L","x: -9.480014<br />y:   0.03612656<br />Identifier_1: NaHCO3 L","x: -9.246623<br />y:  -0.13798188<br />Identifier_1: NaHCO3 L","x: -9.440553<br />y:  -0.26933955<br />Identifier_1: NaHCO3 L","x: -9.359501<br />y:  -0.35181135<br />Identifier_1: NaHCO3 L","x: -9.278409<br />y:  -0.30342202<br />Identifier_1: NaHCO3 L","x:  7.598286<br />y:   5.44654492<br />Identifier_1: NaHCO3 L"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(248,118,109,1)"}},"hoveron":"points","name":"NaHCO3 L","legendgroup":"NaHCO3 L","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[10.327257601598,9.97215116059928,4.92340347275743,9.66309072481126,8.01992583739598,6.16077945628372,5.7833383581412,5.75277042002906,5.35159288402438,5.15844904116612,4.97669175910816,4.68961085273388,4.84909185147737,4.74230506841429,4.67321050643982,7.86862146560608],"y":[4.10471059045198,3.74601880590611,-12.1591117892021,4.08600030014154,3.94000973801211,-12.0716654205766,-12.3953822253477,-12.5745011292579,-12.4737238932092,-12.2835809836757,-11.9091967604582,-11.8485677015563,-12.1712792006585,-11.8112187558857,-11.9692593338478,5.14554270638063],"text":["x: 10.327258<br />y:   4.10471059<br />Identifier_1: NaHCO3 L + NaCl L","x:  9.972151<br />y:   3.74601881<br />Identifier_1: NaHCO3 L + NaCl L","x:  4.923403<br />y: -12.15911179<br />Identifier_1: NaHCO3 L + NaCl L","x:  9.663091<br />y:   4.08600030<br />Identifier_1: NaHCO3 L + NaCl L","x:  8.019926<br />y:   3.94000974<br />Identifier_1: NaHCO3 L + NaCl L","x:  6.160779<br />y: -12.07166542<br />Identifier_1: NaHCO3 L + NaCl L","x:  5.783338<br />y: -12.39538223<br />Identifier_1: NaHCO3 L + NaCl L","x:  5.752770<br />y: -12.57450113<br />Identifier_1: NaHCO3 L + NaCl L","x:  5.351593<br />y: -12.47372389<br />Identifier_1: NaHCO3 L + NaCl L","x:  5.158449<br />y: -12.28358098<br />Identifier_1: NaHCO3 L + NaCl L","x:  4.976692<br />y: -11.90919676<br />Identifier_1: NaHCO3 L + NaCl L","x:  4.689611<br />y: -11.84856770<br />Identifier_1: NaHCO3 L + NaCl L","x:  4.849092<br />y: -12.17127920<br />Identifier_1: NaHCO3 L + NaCl L","x:  4.742305<br />y: -11.81121876<br />Identifier_1: NaHCO3 L + NaCl L","x:  4.673211<br />y: -11.96925933<br />Identifier_1: NaHCO3 L + NaCl L","x:  7.868621<br />y:   5.14554271<br />Identifier_1: NaHCO3 L + NaCl L"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(183,159,0,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(183,159,0,1)"}},"hoveron":"points","name":"NaHCO3 L + NaCl L","legendgroup":"NaHCO3 L + NaCl L","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[9.76349590937271,9.64079213278362,4.98182546464206,9.56813848722826,7.82714141311019,6.68685502381661,6.67168881593628,6.61908808427666,6.4859460520559,6.44617537064551,6.11680476810387,5.98535069657135,5.94130200630185,5.63343317154319,5.30609176621049,7.78283832036656],"y":[3.70127113443431,3.66353175989096,-11.9384162411294,3.98236604106505,3.64165793713852,-11.9533375530156,-11.9911733844026,-12.0153734415358,-12.2854498426474,-12.2548420841026,-12.1226853793825,-12.1931625878698,-12.2711284978253,-12.630133165988,-12.7464315276146,5.15182528065483],"text":["x:  9.763496<br />y:   3.70127113<br />Identifier_1: NaHCO3 L + NaCl U","x:  9.640792<br />y:   3.66353176<br />Identifier_1: NaHCO3 L + NaCl U","x:  4.981825<br />y: -11.93841624<br />Identifier_1: NaHCO3 L + NaCl U","x:  9.568138<br />y:   3.98236604<br />Identifier_1: NaHCO3 L + NaCl U","x:  7.827141<br />y:   3.64165794<br />Identifier_1: NaHCO3 L + NaCl U","x:  6.686855<br />y: -11.95333755<br />Identifier_1: NaHCO3 L + NaCl U","x:  6.671689<br />y: -11.99117338<br />Identifier_1: NaHCO3 L + NaCl U","x:  6.619088<br />y: -12.01537344<br />Identifier_1: NaHCO3 L + NaCl U","x:  6.485946<br />y: -12.28544984<br />Identifier_1: NaHCO3 L + NaCl U","x:  6.446175<br />y: -12.25484208<br />Identifier_1: NaHCO3 L + NaCl U","x:  6.116805<br />y: -12.12268538<br />Identifier_1: NaHCO3 L + NaCl U","x:  5.985351<br />y: -12.19316259<br />Identifier_1: NaHCO3 L + NaCl U","x:  5.941302<br />y: -12.27112850<br />Identifier_1: NaHCO3 L + NaCl U","x:  5.633433<br />y: -12.63013317<br />Identifier_1: NaHCO3 L + NaCl U","x:  5.306092<br />y: -12.74643153<br />Identifier_1: NaHCO3 L + NaCl U","x:  7.782838<br />y:   5.15182528<br />Identifier_1: NaHCO3 L + NaCl U"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,186,56,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,186,56,1)"}},"hoveron":"points","name":"NaHCO3 L + NaCl U","legendgroup":"NaHCO3 L + NaCl U","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[9.01852106462781,9.10675082964558,-9.82117178680755,8.72459836389693,7.43540023211484,-9.83807196535383,-9.76755816113281,-9.66134277682042,-9.81721068093369,-9.89517944025149,-9.84602614505802,-9.68630952073886,-9.59054564768185,-9.43213538746619,-9.26344948406267,7.39454526099014],"y":[3.68900763922402,3.82151139312861,3.86032136563607,3.93632732207876,3.58566077322865,6.67245383292695,6.29744596952841,5.94014975037346,5.21808199488832,4.41824431823177,3.48904880006539,2.94871518591189,2.50868905498803,2.20724560726193,1.68871784304294,5.31898145354475],"text":["x:  9.018521<br />y:   3.68900764<br />Identifier_1: NaHCO3 U","x:  9.106751<br />y:   3.82151139<br />Identifier_1: NaHCO3 U","x: -9.821172<br />y:   3.86032137<br />Identifier_1: NaHCO3 U","x:  8.724598<br />y:   3.93632732<br />Identifier_1: NaHCO3 U","x:  7.435400<br />y:   3.58566077<br />Identifier_1: NaHCO3 U","x: -9.838072<br />y:   6.67245383<br />Identifier_1: NaHCO3 U","x: -9.767558<br />y:   6.29744597<br />Identifier_1: NaHCO3 U","x: -9.661343<br />y:   5.94014975<br />Identifier_1: NaHCO3 U","x: -9.817211<br />y:   5.21808199<br />Identifier_1: NaHCO3 U","x: -9.895179<br />y:   4.41824432<br />Identifier_1: NaHCO3 U","x: -9.846026<br />y:   3.48904880<br />Identifier_1: NaHCO3 U","x: -9.686310<br />y:   2.94871519<br />Identifier_1: NaHCO3 U","x: -9.590546<br />y:   2.50868905<br />Identifier_1: NaHCO3 U","x: -9.432135<br />y:   2.20724561<br />Identifier_1: NaHCO3 U","x: -9.263449<br />y:   1.68871784<br />Identifier_1: NaHCO3 U","x:  7.394545<br />y:   5.31898145<br />Identifier_1: NaHCO3 U"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,191,196,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,191,196,1)"}},"hoveron":"points","name":"NaHCO3 U","legendgroup":"NaHCO3 U","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[9.25498719725089,9.31047434242128,-9.58613038111421,9.12225456781465,7.6693932577299,-9.85515048220016,-9.83531515467042,-9.88122976882916,-9.83601745838686,-9.61292600267066,-9.41269933149145,-9.39803703174131,-9.2574408536526,-9.31093492065312,-9.25163543438447,7.50883592551223],"y":[3.52511481134187,3.72073781343923,2.76748529606576,4.1336774608962,3.57928928271156,6.10802005050461,5.38192135580928,4.59040366895769,3.93625634037222,3.23906313400551,2.62432508257712,2.25750114556471,1.91879293822623,1.67770600583912,1.38884268927744,5.24353214571235],"text":["x:  9.254987<br />y:   3.52511481<br />Identifier_1: NaHCO3 U + NaCl L","x:  9.310474<br />y:   3.72073781<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.586130<br />y:   2.76748530<br />Identifier_1: NaHCO3 U + NaCl L","x:  9.122255<br />y:   4.13367746<br />Identifier_1: NaHCO3 U + NaCl L","x:  7.669393<br />y:   3.57928928<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.855150<br />y:   6.10802005<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.835315<br />y:   5.38192136<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.881230<br />y:   4.59040367<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.836017<br />y:   3.93625634<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.612926<br />y:   3.23906313<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.412699<br />y:   2.62432508<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.398037<br />y:   2.25750115<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.257441<br />y:   1.91879294<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.310935<br />y:   1.67770601<br />Identifier_1: NaHCO3 U + NaCl L","x: -9.251635<br />y:   1.38884269<br />Identifier_1: NaHCO3 U + NaCl L","x:  7.508836<br />y:   5.24353215<br />Identifier_1: NaHCO3 U + NaCl L"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(97,156,255,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(97,156,255,1)"}},"hoveron":"points","name":"NaHCO3 U + NaCl L","legendgroup":"NaHCO3 U + NaCl L","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[7.73099892020039,8.12259460501144,-9.72505750818282,7.91755349625559,7.42933303759377,7.56364824707271,-9.24730132964562,-9.52772783763286,-9.80636210510858,-9.56782521448071,-9.92615095035342,-9.93843609109415,-9.99185679454178,-9.80583620465835,-9.53797517481348,7.37627753172705],"y":[3.8057411765217,3.78095037658213,6.25949810455931,3.5407188973803,3.62044427006491,4.16660158087468,6.8377597976279,6.78211479028233,6.76921607160796,6.46046548154589,5.91289738622867,5.01481917951191,4.16234712094911,3.54997171436972,3.14666338099402,4.96241011783919],"text":["x:  7.730999<br />y:   3.80574118<br />Identifier_1: NaHCO3 U + NaCl U","x:  8.122595<br />y:   3.78095038<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.725058<br />y:   6.25949810<br />Identifier_1: NaHCO3 U + NaCl U","x:  7.917553<br />y:   3.54071890<br />Identifier_1: NaHCO3 U + NaCl U","x:  7.429333<br />y:   3.62044427<br />Identifier_1: NaHCO3 U + NaCl U","x:  7.563648<br />y:   4.16660158<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.247301<br />y:   6.83775980<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.527728<br />y:   6.78211479<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.806362<br />y:   6.76921607<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.567825<br />y:   6.46046548<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.926151<br />y:   5.91289739<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.938436<br />y:   5.01481918<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.991857<br />y:   4.16234712<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.805836<br />y:   3.54997171<br />Identifier_1: NaHCO3 U + NaCl U","x: -9.537975<br />y:   3.14666338<br />Identifier_1: NaHCO3 U + NaCl U","x:  7.376278<br />y:   4.96241012<br />Identifier_1: NaHCO3 U + NaCl U"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(245,100,227,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(245,100,227,1)"}},"hoveron":"points","name":"NaHCO3 U + NaCl U","legendgroup":"NaHCO3 U + NaCl U","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":43.7625570776256,"r":7.30593607305936,"b":40.1826484018265,"l":43.1050228310502},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"title":{"text":"UMAP of NaHCO3 + NaCl (not quality-checked)","font":{"color":"rgba(0,0,0,1)","family":"","size":17.5342465753425},"x":0,"xref":"paper"},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-11.0127554751478,11.4470154981846],"tickmode":"array","ticktext":["-10","-5","0","5","10"],"tickvals":[-10,-5,0,5,10],"categoryorder":"array","categoryarray":["-10","-5","0","5","10"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":{"text":"x","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-13.7256410938767,7.81696936389003],"tickmode":"array","ticktext":["-10","-5","0","5"],"tickvals":[-10,-5,0,5],"categoryorder":"array","categoryarray":["-10","-5","0","5"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":{"text":"y","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"y":0.93503937007874},"annotations":[{"text":"Identifier_1","x":1.02,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"left","yanchor":"bottom","legendTitle":true}],"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"11ef640450a3e":{"x":{},"y":{},"colour":{},"type":"scatter"}},"cur_data":"11ef640450a3e","visdat":{"11ef640450a3e":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```
<br><br>


## Quality-controlled data
UMAP for different reactions that have been quality-checked.<br>

EXAMPLE of interactive plot with multiple variables


```{.r .details .show summary='Source'}
library(gapminder) 
p <- gapminder %>%
  filter(year==1977) %>%
  ggplot( aes(gdpPercap, lifeExp, size = pop, color=continent)) +
  geom_point() +
  scale_x_log10() +
  theme_bw()
 
ggplotly(p)
```

```{=html}
<div id="htmlwidget-72d4cdf0f6a211518560" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-72d4cdf0f6a211518560">{"x":{"data":[{"x":[3.69111835304969,3.47837128685849,3.0124834262036,3.50716177034781,2.87121498358872,2.74515544499733,3.25125676748187,3.04507811557752,3.05460729073385,3.06915101836112,2.90078062155373,3.51310821042633,3.40101028408395,3.48879895792208,3.44490216214739,2.98162238872728,2.7039391611549,2.74570576546197,4.3373708614867,2.94682314869669,2.99704718642708,2.94185210799775,2.88350583493768,3.10298675433845,2.87237164152654,2.8063987207202,4.34145849934605,3.18871158763867,2.82166002181421,2.83657428165302,3.17536457562432,3.56948896052176,3.37486193965761,2.70098023955392,3.58843821514052,2.90789326381397,3.29709308978461,3.6354640490023,2.82612704524949,3.23994022438144,3.19361683040004,3.12978175413031,3.16166517152032,3.90464260358686,3.34301221488094,3.57765383934322,2.9833972607425,3.185478974435,3.49427662628827,2.92620510614229,3.20103869689067,2.83606300605042],"y":[58.014,39.483,49.19,59.319,46.137,45.91,49.355,46.775,47.383,50.939,47.804,55.625,52.374,46.519,53.319,42.024,44.535,44.51,52.79,41.842,51.756,40.762,37.465,56.155,52.208,43.764,57.442,46.881,43.767,41.714,50.852,64.93,55.73,42.495,56.437,41.291,44.514,67.064,45,58.55,48.879,36.788,41.974,55.527,47.8,52.537,49.919,52.887,59.837,50.35,51.386,57.674],"text":["gdpPercap:  4910.4168<br />lifeExp: 58.01400<br />pop:  17152804<br />continent: Africa","gdpPercap:  3008.6474<br />lifeExp: 39.48300<br />pop:   6162675<br />continent: Africa","gdpPercap:  1029.1613<br />lifeExp: 49.19000<br />pop:   3168267<br />continent: Africa","gdpPercap:  3214.8578<br />lifeExp: 59.31900<br />pop:    781472<br />continent: Africa","gdpPercap:   743.3870<br />lifeExp: 46.13700<br />pop:   5889574<br />continent: Africa","gdpPercap:   556.1033<br />lifeExp: 45.91000<br />pop:   3834415<br />continent: Africa","gdpPercap:  1783.4329<br />lifeExp: 49.35500<br />pop:   7959865<br />continent: Africa","gdpPercap:  1109.3743<br />lifeExp: 46.77500<br />pop:   2167533<br />continent: Africa","gdpPercap:  1133.9850<br />lifeExp: 47.38300<br />pop:   4388260<br />continent: Africa","gdpPercap:  1172.6030<br />lifeExp: 50.93900<br />pop:    304739<br />continent: Africa","gdpPercap:   795.7573<br />lifeExp: 47.80400<br />pop:  26480870<br />continent: Africa","gdpPercap:  3259.1790<br />lifeExp: 55.62500<br />pop:   1536769<br />continent: Africa","gdpPercap:  2517.7365<br />lifeExp: 52.37400<br />pop:   7459574<br />continent: Africa","gdpPercap:  3081.7610<br />lifeExp: 46.51900<br />pop:    228694<br />continent: Africa","gdpPercap:  2785.4936<br />lifeExp: 53.31900<br />pop:  38783863<br />continent: Africa","gdpPercap:   958.5668<br />lifeExp: 42.02400<br />pop:    192675<br />continent: Africa","gdpPercap:   505.7538<br />lifeExp: 44.53500<br />pop:   2512642<br />continent: Africa","gdpPercap:   556.8084<br />lifeExp: 44.51000<br />pop:  34617799<br />continent: Africa","gdpPercap: 21745.5733<br />lifeExp: 52.79000<br />pop:    706367<br />continent: Africa","gdpPercap:   884.7553<br />lifeExp: 41.84200<br />pop:    608274<br />continent: Africa","gdpPercap:   993.2240<br />lifeExp: 51.75600<br />pop:  10538093<br />continent: Africa","gdpPercap:   874.6859<br />lifeExp: 40.76200<br />pop:   4227026<br />continent: Africa","gdpPercap:   764.7260<br />lifeExp: 37.46500<br />pop:    745228<br />continent: Africa","gdpPercap:  1267.6132<br />lifeExp: 56.15500<br />pop:  14500404<br />continent: Africa","gdpPercap:   745.3695<br />lifeExp: 52.20800<br />pop:   1251524<br />continent: Africa","gdpPercap:   640.3224<br />lifeExp: 43.76400<br />pop:   1703617<br />continent: Africa","gdpPercap: 21951.2118<br />lifeExp: 57.44200<br />pop:   2721783<br />continent: Africa","gdpPercap:  1544.2286<br />lifeExp: 46.88100<br />pop:   8007166<br />continent: Africa","gdpPercap:   663.2237<br />lifeExp: 43.76700<br />pop:   5637246<br />continent: Africa","gdpPercap:   686.3953<br />lifeExp: 41.71400<br />pop:   6491649<br />continent: Africa","gdpPercap:  1497.4922<br />lifeExp: 50.85200<br />pop:   1456688<br />continent: Africa","gdpPercap:  3710.9830<br />lifeExp: 64.93000<br />pop:    913025<br />continent: Africa","gdpPercap:  2370.6200<br />lifeExp: 55.73000<br />pop:  18396941<br />continent: Africa","gdpPercap:   502.3197<br />lifeExp: 42.49500<br />pop:  11127868<br />continent: Africa","gdpPercap:  3876.4860<br />lifeExp: 56.43700<br />pop:    977026<br />continent: Africa","gdpPercap:   808.8971<br />lifeExp: 41.29100<br />pop:   5682086<br />continent: Africa","gdpPercap:  1981.9518<br />lifeExp: 44.51400<br />pop:  62209173<br />continent: Africa","gdpPercap:  4319.8041<br />lifeExp: 67.06400<br />pop:    492095<br />continent: Africa","gdpPercap:   670.0806<br />lifeExp: 45.00000<br />pop:   4657072<br />continent: Africa","gdpPercap:  1737.5617<br />lifeExp: 58.55000<br />pop:     86796<br />continent: Africa","gdpPercap:  1561.7691<br />lifeExp: 48.87900<br />pop:   5260855<br />continent: Africa","gdpPercap:  1348.2852<br />lifeExp: 36.78800<br />pop:   3140897<br />continent: Africa","gdpPercap:  1450.9925<br />lifeExp: 41.97400<br />pop:   4353666<br />continent: Africa","gdpPercap:  8028.6514<br />lifeExp: 55.52700<br />pop:  27129932<br />continent: Africa","gdpPercap:  2202.9884<br />lifeExp: 47.80000<br />pop:  17104986<br />continent: Africa","gdpPercap:  3781.4106<br />lifeExp: 52.53700<br />pop:    551425<br />continent: Africa","gdpPercap:   962.4923<br />lifeExp: 49.91900<br />pop:  17129565<br />continent: Africa","gdpPercap:  1532.7770<br />lifeExp: 52.88700<br />pop:   2308582<br />continent: Africa","gdpPercap:  3120.8768<br />lifeExp: 59.83700<br />pop:   6005061<br />continent: Africa","gdpPercap:   843.7331<br />lifeExp: 50.35000<br />pop:  11457758<br />continent: Africa","gdpPercap:  1588.6883<br />lifeExp: 51.38600<br />pop:   5216550<br />continent: Africa","gdpPercap:   685.5877<br />lifeExp: 57.67400<br />pop:   6642107<br />continent: Africa"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","opacity":1,"size":[6.32127781530702,5.29612899121632,4.85958274371316,4.29233940761101,5.26165272498491,4.97061767194295,5.5059176412168,4.66704214356912,5.05559917855771,4.06676322554192,6.94049569932146,4.52040583610645,5.45016604603988,4.01129639562885,7.60694585246374,3.97973112975184,4.73782100329548,7.39505365224446,4.26382530133105,4.22383634289638,5.76860640347983,5.03145493864573,4.27878249624843,6.11542223529683,4.44354450304809,4.56187177172394,4.77827606357594,5.51109590728372,5.22907021838905,5.3366454219231,4.49965617393791,4.33879188894232,6.41229679974379,5.82395887424691,4.36004866819238,5.23491359575959,8.62896025232351,4.17122854948148,5.09486783222961,3.77952755905512,5.17905878226218,4.854775454452,5.05045749460444,6.97912565512447,6.31771439639559,4.19891947998289,6.31954665895726,4.69663035294985,5.27632870617264,5.85427629024999,5.17305387000117,5.35482857766587],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(248,118,109,1)"}},"hoveron":"points","name":"Africa","legendgroup":"Africa","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[4.00341859740172,3.54999558616307,3.82348196645535,4.34421307668719,3.67731159008544,3.58158649939174,3.77282591203612,3.80485437036137,3.42845697609622,3.82475196835622,3.71087205749874,3.68841917661215,3.27283885748582,3.50558514939778,3.82283441748743,3.88507437266432,3.7392851792172,3.72850897552038,3.51166593362054,3.7980489038637,3.98991789679654,3.89760258369119,4.38152357914239,3.81320321246222,4.11872592974958],"y":[68.481,50.023,61.489,74.21,67.052,63.837,70.75,72.649,61.788,61.31,56.696,56.029,49.923,57.402,70.11,65.032,57.47,68.681,66.353,58.447,73.44,68.3,73.38,69.481,67.456],"text":["gdpPercap: 10079.0267<br />lifeExp: 68.48100<br />pop:  26983828<br />continent: Americas","gdpPercap:  3548.0978<br />lifeExp: 50.02300<br />pop:   5079716<br />continent: Americas","gdpPercap:  6660.1187<br />lifeExp: 61.48900<br />pop: 114313951<br />continent: Americas","gdpPercap: 22090.8831<br />lifeExp: 74.21000<br />pop:  23796400<br />continent: Americas","gdpPercap:  4756.7638<br />lifeExp: 67.05200<br />pop:  10599793<br />continent: Americas","gdpPercap:  3815.8079<br />lifeExp: 63.83700<br />pop:  25094412<br />continent: Americas","gdpPercap:  5926.8770<br />lifeExp: 70.75000<br />pop:   2108457<br />continent: Americas","gdpPercap:  6380.4950<br />lifeExp: 72.64900<br />pop:   9537988<br />continent: Americas","gdpPercap:  2681.9889<br />lifeExp: 61.78800<br />pop:   5302800<br />continent: Americas","gdpPercap:  6679.6233<br />lifeExp: 61.31000<br />pop:   7278866<br />continent: Americas","gdpPercap:  5138.9224<br />lifeExp: 56.69600<br />pop:   4282586<br />continent: Americas","gdpPercap:  4879.9927<br />lifeExp: 56.02900<br />pop:   5703430<br />continent: Americas","gdpPercap:  1874.2989<br />lifeExp: 49.92300<br />pop:   4908554<br />continent: Americas","gdpPercap:  3203.2081<br />lifeExp: 57.40200<br />pop:   3055235<br />continent: Americas","gdpPercap:  6650.1956<br />lifeExp: 70.11000<br />pop:   2156814<br />continent: Americas","gdpPercap:  7674.9291<br />lifeExp: 65.03200<br />pop:  63759976<br />continent: Americas","gdpPercap:  5486.3711<br />lifeExp: 57.47000<br />pop:   2554598<br />continent: Americas","gdpPercap:  5351.9121<br />lifeExp: 68.68100<br />pop:   1839782<br />continent: Americas","gdpPercap:  3248.3733<br />lifeExp: 66.35300<br />pop:   2984494<br />continent: Americas","gdpPercap:  6281.2909<br />lifeExp: 58.44700<br />pop:  15990099<br />continent: Americas","gdpPercap:  9770.5249<br />lifeExp: 73.44000<br />pop:   3080828<br />continent: Americas","gdpPercap:  7899.5542<br />lifeExp: 68.30000<br />pop:   1039009<br />continent: Americas","gdpPercap: 24072.6321<br />lifeExp: 73.38000<br />pop: 220239000<br />continent: Americas","gdpPercap:  6504.3397<br />lifeExp: 69.48100<br />pop:   2873520<br />continent: Americas","gdpPercap: 13143.9510<br />lifeExp: 67.45600<br />pop:  13503563<br />continent: Americas"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(163,165,0,1)","opacity":1,"size":[6.97047083068522,5.15434238640086,10.3553727984455,6.77543956769576,5.77446910007444,6.85635435789675,4.65435232683653,5.6710443310913,5.18472018302092,5.42956521861593,5.0398270921093,5.23768683822733,5.13057191783969,4.83958883002819,4.66475315834704,8.68911697083458,4.74607251669638,4.59414954550014,4.82688148286866,6.23316601826461,4.84414878929001,4.37991825947829,12.9086371873162,4.80663028657321,6.03320025163606],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(163,165,0,1)"}},"hoveron":"points","name":"Americas","legendgroup":"Americas","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[2.89548517717506,4.28645875933068,2.81946314412135,2.72013629197205,2.86995736499049,4.048680298964,2.91027070197067,3.14072860861607,4.07513053543838,4.16696961434647,4.12406772903255,4.22037948985684,3.45520305365965,3.61345080718246,3.66812684894598,4.7728017852095,3.93750268823479,3.58296303108212,3.21682849824462,2.56937390961505,2.84142982784579,4.07365765189317,3.07037821749517,3.37533512418369,4.53361654092087,4.04960907918865,3.12993971717703,3.50453672644772,3.74791804661293,3.29252733979409,2.85341657088882,3.56618184839601,3.26239535810376],"y":[38.438,65.593,46.923,31.22,63.96736,73.6,54.208,52.702,57.702,60.413,73.06,75.38,61.134,67.159,64.766,69.343,66.099,65.256,55.491,56.059,46.748,57.367,54.043,60.06,58.69,70.795,65.949,61.195,70.59,62.494,55.764,60.765,44.175],"text":["gdpPercap:   786.1134<br />lifeExp: 38.43800<br />pop:  14880372<br />continent: Asia","gdpPercap: 19340.1020<br />lifeExp: 65.59300<br />pop:    297410<br />continent: Asia","gdpPercap:   659.8772<br />lifeExp: 46.92300<br />pop:  80428306<br />continent: Asia","gdpPercap:   524.9722<br />lifeExp: 31.22000<br />pop:   6978607<br />continent: Asia","gdpPercap:   741.2375<br />lifeExp: 63.96736<br />pop: 943455000<br />continent: Asia","gdpPercap: 11186.1413<br />lifeExp: 73.60000<br />pop:   4583700<br />continent: Asia","gdpPercap:   813.3373<br />lifeExp: 54.20800<br />pop: 634000000<br />continent: Asia","gdpPercap:  1382.7021<br />lifeExp: 52.70200<br />pop: 136725000<br />continent: Asia","gdpPercap: 11888.5951<br />lifeExp: 57.70200<br />pop:  35480679<br />continent: Asia","gdpPercap: 14688.2351<br />lifeExp: 60.41300<br />pop:  11882916<br />continent: Asia","gdpPercap: 13306.6192<br />lifeExp: 73.06000<br />pop:   3495918<br />continent: Asia","gdpPercap: 16610.3770<br />lifeExp: 75.38000<br />pop: 113872473<br />continent: Asia","gdpPercap:  2852.3516<br />lifeExp: 61.13400<br />pop:   1937652<br />continent: Asia","gdpPercap:  4106.3012<br />lifeExp: 67.15900<br />pop:  16325320<br />continent: Asia","gdpPercap:  4657.2210<br />lifeExp: 64.76600<br />pop:  36436000<br />continent: Asia","gdpPercap: 59265.4771<br />lifeExp: 69.34300<br />pop:   1140357<br />continent: Asia","gdpPercap:  8659.6968<br />lifeExp: 66.09900<br />pop:   3115787<br />continent: Asia","gdpPercap:  3827.9216<br />lifeExp: 65.25600<br />pop:  12845381<br />continent: Asia","gdpPercap:  1647.5117<br />lifeExp: 55.49100<br />pop:   1528000<br />continent: Asia","gdpPercap:   371.0000<br />lifeExp: 56.05900<br />pop:  31528087<br />continent: Asia","gdpPercap:   694.1124<br />lifeExp: 46.74800<br />pop:  13933198<br />continent: Asia","gdpPercap: 11848.3439<br />lifeExp: 57.36700<br />pop:   1004533<br />continent: Asia","gdpPercap:  1175.9212<br />lifeExp: 54.04300<br />pop:  78152686<br />continent: Asia","gdpPercap:  2373.2043<br />lifeExp: 60.06000<br />pop:  46850962<br />continent: Asia","gdpPercap: 34167.7626<br />lifeExp: 58.69000<br />pop:   8128505<br />continent: Asia","gdpPercap: 11210.0895<br />lifeExp: 70.79500<br />pop:   2325300<br />continent: Asia","gdpPercap:  1348.7757<br />lifeExp: 65.94900<br />pop:  14116836<br />continent: Asia","gdpPercap:  3195.4846<br />lifeExp: 61.19500<br />pop:   7932503<br />continent: Asia","gdpPercap:  5596.5198<br />lifeExp: 70.59000<br />pop:  16785196<br />continent: Asia","gdpPercap:  1961.2246<br />lifeExp: 62.49400<br />pop:  44148285<br />continent: Asia","gdpPercap:   713.5371<br />lifeExp: 55.76400<br />pop:  50533506<br />continent: Asia","gdpPercap:  3682.8315<br />lifeExp: 60.76500<br />pop:   1261091<br />continent: Asia","gdpPercap:  1829.7652<br />lifeExp: 44.17500<br />pop:   8403990<br />continent: Asia"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,191,125,1)","opacity":1,"size":[6.14601109401926,4.06189233680308,9.29441811927379,5.39475462377927,22.6771653543307,5.08426676230058,19.2706113323364,10.9715821633898,7.43994831281003,5.89270757284543,4.9155533082844,10.342652961876,4.61658106384004,6.25889092142746,7.4890188516279,4.41106157328659,4.85034613145665,5.97722629927717,4.51816213431315,7.2295116202132,6.06899974449858,4.36894912815136,9.21575428359059,7.98702151037774,5.52430918996644,4.70007429187427,6.08413180764231,5.50291508915639,6.29375367739212,7.86362869673143,8.14954658274176,4.44626602186696,5.5539430462659],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,191,125,1)"}},"hoveron":"points","name":"Asia","legendgroup":"Asia","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[3.54814411807184,4.29555439638822,4.28144187754483,3.54758782074795,3.88151249733011,4.05328536317461,4.17026642863828,4.31011744274222,4.19327554039949,4.26227627214274,4.31202751232983,4.15215143685527,4.06725083957618,4.29347221889702,4.04731308089082,4.15399722211687,3.98208706693409,4.32652140433408,4.36756741361793,3.97809563416623,4.00742708886894,3.97110865215781,4.11329709452241,4.03832857580155,4.18403893629866,4.12178698250771,4.27544324046184,4.43107881405245,3.63033859911073,4.24126620200209],"y":[68.93,72.17,72.8,69.86,70.81,70.64,70.71,74.69,72.52,73.83,72.5,73.68,69.95,76.11,72.03,73.48,73.066,75.24,75.37,70.67,70.41,69.46,70.3,70.45,70.97,74.39,75.44,75.39,59.507,72.76],"text":["gdpPercap:  3533.0039<br />lifeExp: 68.93000<br />pop:   2509048<br />continent: Europe","gdpPercap: 19749.4223<br />lifeExp: 72.17000<br />pop:   7568430<br />continent: Europe","gdpPercap: 19117.9745<br />lifeExp: 72.80000<br />pop:   9821800<br />continent: Europe","gdpPercap:  3528.4813<br />lifeExp: 69.86000<br />pop:   4086000<br />continent: Europe","gdpPercap:  7612.2404<br />lifeExp: 70.81000<br />pop:   8797022<br />continent: Europe","gdpPercap: 11305.3852<br />lifeExp: 70.64000<br />pop:   4318673<br />continent: Europe","gdpPercap: 14800.1606<br />lifeExp: 70.71000<br />pop:  10161915<br />continent: Europe","gdpPercap: 20422.9015<br />lifeExp: 74.69000<br />pop:   5088419<br />continent: Europe","gdpPercap: 15605.4228<br />lifeExp: 72.52000<br />pop:   4738902<br />continent: Europe","gdpPercap: 18292.6351<br />lifeExp: 73.83000<br />pop:  53165019<br />continent: Europe","gdpPercap: 20512.9212<br />lifeExp: 72.50000<br />pop:  78160773<br />continent: Europe","gdpPercap: 14195.5243<br />lifeExp: 73.68000<br />pop:   9308479<br />continent: Europe","gdpPercap: 11674.8374<br />lifeExp: 69.95000<br />pop:  10637171<br />continent: Europe","gdpPercap: 19654.9625<br />lifeExp: 76.11000<br />pop:    221823<br />continent: Europe","gdpPercap: 11150.9811<br />lifeExp: 72.03000<br />pop:   3271900<br />continent: Europe","gdpPercap: 14255.9847<br />lifeExp: 73.48000<br />pop:  56059245<br />continent: Europe","gdpPercap:  9595.9299<br />lifeExp: 73.06600<br />pop:    560073<br />continent: Europe","gdpPercap: 21209.0592<br />lifeExp: 75.24000<br />pop:  13852989<br />continent: Europe","gdpPercap: 23311.3494<br />lifeExp: 75.37000<br />pop:   4043205<br />continent: Europe","gdpPercap:  9508.1415<br />lifeExp: 70.67000<br />pop:  34621254<br />continent: Europe","gdpPercap: 10172.4857<br />lifeExp: 70.41000<br />pop:   9662600<br />continent: Europe","gdpPercap:  9356.3972<br />lifeExp: 69.46000<br />pop:  21658597<br />continent: Europe","gdpPercap: 12980.6696<br />lifeExp: 70.30000<br />pop:   8686367<br />continent: Europe","gdpPercap: 10922.6640<br />lifeExp: 70.45000<br />pop:   4827803<br />continent: Europe","gdpPercap: 15277.0302<br />lifeExp: 70.97000<br />pop:   1746919<br />continent: Europe","gdpPercap: 13236.9212<br />lifeExp: 74.39000<br />pop:  36439000<br />continent: Europe","gdpPercap: 18855.7252<br />lifeExp: 75.44000<br />pop:   8251648<br />continent: Europe","gdpPercap: 26982.2905<br />lifeExp: 75.39000<br />pop:   6316424<br />continent: Europe","gdpPercap:  4269.1223<br />lifeExp: 59.50700<br />pop:  42404033<br />continent: Europe","gdpPercap: 17428.7485<br />lifeExp: 72.76000<br />pop:  56179000<br />continent: Europe"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,176,246,1)","opacity":1,"size":[4.73711086274366,5.46245399673978,5.69923465730081,5.00994841135964,5.59538443706805,5.04523525790558,5.73248147708556,5.15554006271282,5.10659105931128,8.26207702648979,9.21603585103979,5.64793676344929,5.77801236958701,4.00561540333621,4.87759423629259,8.38266622872152,4.20280449359801,6.06235892227666,5.00334741211677,7.39523452347783,5.68347312501541,6.63718460047671,5.58381321630558,5.11921103857739,4.57227899402677,7.48917192574021,5.53761739047554,5.31519777193933,7.78197424302696,8.38758788519898],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,176,246,1)"}},"hoveron":"points","name":"Europe","legendgroup":"Europe","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[4.26326190559163,4.21041798943211],"y":[73.49,72.22],"text":["gdpPercap: 18334.1975<br />lifeExp: 73.49000<br />pop:  14074100<br />continent: Oceania","gdpPercap: 16233.7177<br />lifeExp: 72.22000<br />pop:   3164900<br />continent: Oceania"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(231,107,243,1)","opacity":1,"size":[6.08061917751821,4.85899251590851],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(231,107,243,1)"}},"hoveron":"points","name":"Oceania","legendgroup":"Oceania","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.2283105022831,"r":7.30593607305936,"b":40.1826484018265,"l":37.2602739726027},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[2.45920251583532,4.88297317898922],"tickmode":"array","ticktext":["300","1000","3000","10000","30000"],"tickvals":[2.47712125471966,3,3.47712125471966,4,4.47712125471966],"categoryorder":"array","categoryarray":["300","1000","3000","10000","30000"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":{"text":"gdpPercap","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[28.9755,78.3545],"tickmode":"array","ticktext":["30","40","50","60","70"],"tickvals":[30,40,50,60,70],"categoryorder":"array","categoryarray":["30","40","50","60","70"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":{"text":"lifeExp","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"y":0.87007874015748},"annotations":[{"text":"continent<br />pop","x":1.02,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"left","yanchor":"bottom","legendTitle":true}],"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"11ef63bf5bdd1":{"x":{},"y":{},"size":{},"colour":{},"type":"scatter"}},"cur_data":"11ef63bf5bdd1","visdat":{"11ef63bf5bdd1":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```
<br>


# References


