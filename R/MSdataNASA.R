setwd("~lily/Desktop/MASPEX")
#install.packages("isoreader")
library(isoreader)

# NASA example data
# read continuous flow
msdat<-iso_read_continuous_flow("170525_NaHCO3 L + NaCl U_.dxf")
class(msdat)
# [1] "continuous_flow" "iso_file"

# file summary
#msdat %>% iso_get_data_summary() %>% knitr::kable()
summary(msdat)
#   Length Class           Mode
# version            1     package_version list
# read_options       4     -none-          list
# file_info         16     tbl_df          list
# method_info        3     -none-          list
# raw_data           5     tbl_df          list
# vendor_data_table 39     tbl_df          list

read_summ<-function(filename){
  msdat<-iso_read_continuous_flow(filename)
  summ<-summary(msdat)
  print(summ)
}

summ<-read_summ("170525_NaHCO3 L + NaCl U_.dxf")
class(summ)
str(msdat)

msdat$version
# [1] ‘1.3.0’


msdat$method_info
# $standards
# A tibble: 2 x 5
#standard gas   delta_name delta_value reference
#<chr>    <chr> <chr>            <dbl> <chr>
#  1 CO2_zero CO2   d 13C/12C        -36.9 VPDB
#. 2 CO2_zero CO2   d 18O/16O        -40   VSMOW
#
#$reference_ratios
# A tibble: 6 x 4
#reference element ratio_name ratio_value
#<chr>     <chr>   <chr>            <dbl>
# 1 VPDB      C       R 13C/12C     0.0112
# 2 VPDB      O       R 18O/16O     0.00207
# 3 VPDB      O       R 17O/16O     0.000386
# 4 VSMOW     H       R 2H/1H       0.000156
# 5 VSMOW     O       R 17O/16O     0.000380
# 6 VSMOW     O       R 18O/16O     0.00201
#
#$resistors
# A tibble: 3 x 3
#cup        R.Ohm mass
#<int>        <dbl> <chr>
# 1     1    300000000 44
# 2     2  30000000000 45
# 3     3 100000000000 46


msdat$read_options
# $file_info
#[1] TRUE
#
#$method_info
#[1] TRUE
#
#$raw_data
#[1] TRUE
#
#$vendor_data_table
#[1] TRUE


msdat$file_info
# # A tibble: 1 x 16
#file_id      file_root file_path      file_subpath file_datetime       file_size Row   `Peak Center` `AS Sample` `AS Method`
#<chr>        <chr>     <chr>          <chr>        <dttm>                  <int> <lis> <list>        <list>      <list>
#  1 170525_NaHC… .         170525_NaHCO3… NA           2017-05-25 21:32:18    436817 <chr… <chr [1]>     <chr [1]>   <chr [1]>
  # … with 6 more variables: Identifier 1 <list>, Analysis <list>, Preparation <list>, Method <list>, measurement_info <list>,
  #   MS_integration_time.s <list>

# can access by index
msdat$file_info[16]


msdat$raw_data
# # A tibble: 4,297 x 5
#       tp  time.s. v44.mV v45.mV v46.mV
#.    <int>  <dbl>  <dbl>  <dbl>  <dbl>
# 1     1  0.209   1.22  0.578   2.51
# 2     2  0.418   1.21  0.633   2.64
# 3     3  0.627   1.21  0.463   2.26
# 4     4  0.836   1.21  0.593   2.53
# 5     5  1.04    1.21  0.559   2.55
# 6     6  1.25    1.22  0.572   2.53
# 7     7  1.46    1.22  0.656   2.59
# 8     8  1.67    1.22  0.860   3.25
# 9     9  1.88    1.22  0.461   2.26
# 10    10  2.09    1.21  0.522   2.23
# … with 4,287 more rows

msdat$vendor_data_table
# # A tibble: 16 x 39
#      Nr.    Start       Rt      End `Ampl 44` `Ampl 45` `Ampl 46`  `BGD 44`  `BGD 45` `BGD 46` `rIntensity 44` `rIntensity 45`
#<.   int> <dbl[s]> <dbl[s]> <dbl[s]> <dbl[mV]> <dbl[mV]> <dbl[mV]> <dbl[mV]> <dbl[mV]> <dbl[mV>      <dbl[mVs]>      <dbl[mVs]>
# 1     1   27.170   47.443   50.369 2061.8995 2372.8830 2811.8782  1.296127 0.6860018 2.513521       40127.019       46179.141
#2     2   67.089   87.153   90.079 2062.2813 2373.3475 2814.5886  1.496863 0.8692536 2.886713       40120.612       46171.250
#3     3  131.879  134.178  138.358  692.0574  816.8932  953.1559  1.512164 0.9017184 2.639453        1070.433        1262.219
#4     4  166.573  186.846  189.563 2062.3738 2372.9537 2813.6276  1.385960 0.8387023 2.999882       40106.546       46154.659
#5     5  206.492  226.556  229.482 2059.6704 2370.0244 2809.4172  1.519815 0.8387023 2.756344       39909.064       45925.618
#6     6  271.282  275.253  281.941 1119.1845 1318.9636 1583.0331  1.506426 0.8024275 2.940414        4341.673        5117.680
#7     7  331.056  334.818  341.506 1064.7337 1255.2545 1505.7065  1.569550 1.1252670 3.396016        4130.948        4869.138
#8     8  390.621  394.592  401.280 1014.9099 1196.0384 1434.5239  1.460528 0.8635249 2.846445        3935.767        4638.832
#9     9  450.395  454.366  460.845  967.3212 1139.9215 1368.5723  1.472001 1.0048699 3.090960        3750.687        4420.134
#10    10  510.169  513.931  520.410  922.2640 1086.8234 1304.9727  1.512164 0.8520680 2.670107        3575.083        4213.020
#11    11  569.943  573.705  580.184  881.1731 1038.1846 1245.1615  1.575290 1.0813058 3.133147        3411.304        4019.923
#12    12  629.508  633.479  639.749  842.2562  992.4550 1191.2397  1.410814 0.7127176 2.606886        3259.319        3840.505
#13    13  689.282  693.253  699.314  805.8029  949.4437 1140.2887  1.410814 0.8730728 3.059364        3115.831        3670.764
#14    14  749.056  752.818  759.088  771.7275  909.0003 1091.1399  1.470089 0.8406116 2.842610        2977.501        3507.271
#15    15  808.621  812.592  818.653  739.3136  870.7118 1045.6774  1.357287 0.7795198 2.767844        2854.095        3361.512
#16    16  843.315  863.588  866.305 2056.0245 2366.2184 2805.1346  1.428021 0.7604316 2.647116       39958.264       45985.169
# … with 27 more variables: rIntensity 46 <dbl[mVs]>, rIntensity All <dbl[mVs]>, Intensity 44 <dbl[Vs]>,
#   Intensity 45 <dbl[Vs]>, Intensity 46 <dbl[Vs]>, Intensity All <dbl[Vs]>, List First Peak <int>, rR 45CO2/44CO2 <dbl>,
#   rR 46CO2/44CO2 <dbl>, Is Ref.? <int>, R 45CO2/44CO2 <dbl>, Ref. Name <chr>, rd 45CO2/44CO2 <dbl[permil]>,
#   d 45CO2/44CO2 <dbl[permil]>, R 46CO2/44CO2 <dbl>, rd 46CO2/44CO2 <dbl[permil]>, d 46CO2/44CO2 <dbl[permil]>,
#   R 13C/12C <dbl>, d 13C/12C <dbl[permil]>, AT% 13C/12C <dbl[%]>, R 18O/16O <dbl>, d 18O/16O <dbl[permil]>,
#   AT% 18O/16O <dbl[%]>, R 17O/16O <dbl>, d 17O/16O <dbl>, Rps 45CO2/44CO2 <dbl>, Rps 46CO2/44CO2 <dbl>





# if trouble reading any of the .dxf files
#msdat %>% iso_get_problems_summary() %>% knitr::kable()
#msdat %>% iso_get_problems() %>% knitr::kable()
iso_get_problems_summary(msdat) # no problems
iso_get_problems(msdat) # no problems


# all file information
#allInfo<-msdat %>% iso_get_file_info(select = c(-file_root)) %>% knitr::kable()
allInfo<-iso_get_file_info(msdat)
allInfo
# Info: aggregating file info from 1 data file(s)
# A tibble: 1 x 16
#file_id      file_root file_path      file_subpath file_datetime       file_size Row   `Peak Center` `AS Sample` `AS Method`
#<chr>        <chr>     <chr>          <chr>        <dttm>                  <int> <chr> <chr>         <chr>       <chr>
#  1 170525_NaHC… .         170525_NaHCO3… NA           2017-05-25 21:32:18    436817 8     1             8           >Internal …
# … with 6 more variables: Identifier 1 <chr>, Analysis <chr>, Preparation <chr>, Method <chr>, measurement_info <list>,
#   MS_integration_time.s <dbl>

(allInfo$measurement_info)

# select file information
msdat %>%
  iso_get_file_info(
    select = c(
      # rename sample id columns from the different file types to a new ID column
      ID = `Identifier 1`, ID = `Name`,
      # select columns without renaming
      Analysis, `Peak Center`, `H3 Factor`,
      # select the time stamp and rename it to `Date & Time`
      `Date & Time` = file_datetime
    ),
    # explicitly allow for file specific rename (for the new ID column)
    file_specific = TRUE
  ) %>% knitr::kable()


# Resistors
# some IRMS data files contain resistor information useful for downstream calcs (i.e. signal conversion)
#msdat %>% iso_get_resistors() %>% knitr::kable()
iso_get_resistors(msdat)
#Info: aggregating resistors info from 1 data file(s)
# A tibble: 3 x 4
#file_id                           cup      R.Ohm mass
#<chr>                            <int>     <dbl> <chr>
# 1 170525_NaHCO3 L + NaCl U_.dxf   1    300000000 44
#2 170525_NaHCO3 L + NaCl U_.dxf    2  30000000000 45
#3 170525_NaHCO3 L + NaCl U_.dxf    3 100000000000 46


# Reference values
# isotopic reference values for different gases
# reference delta values without ratio values
msdat %>% iso_get_standards(file_id:reference) %>% knitr::kable()
# reference values with ratios
msdat %>% iso_get_standards() %>% knitr::kable()


# get raw data with default selections (all raw data, no additional file info)
#msdat %>% iso_get_raw_data() %>% head(n=10) %>% knitr::kable()
rawdat<-iso_get_raw_data(msdat)
rawdat
# # A tibble: 4,297 x 6
#file_id                          tp time.s v44.mV v45.mV v46.mV
#<chr>                         <int>  <dbl>  <dbl>  <dbl>  <dbl>
# 1 170525_NaHCO3 L + NaCl U_.dxf     1  0.209   1.22  0.578   2.51
# 2 170525_NaHCO3 L + NaCl U_.dxf     2  0.418   1.21  0.633   2.64
# 3 170525_NaHCO3 L + NaCl U_.dxf     3  0.627   1.21  0.463   2.26
# 4 170525_NaHCO3 L + NaCl U_.dxf     4  0.836   1.21  0.593   2.53
# 5 170525_NaHCO3 L + NaCl U_.dxf     5  1.04    1.21  0.559   2.55
# 6 170525_NaHCO3 L + NaCl U_.dxf     6  1.25    1.22  0.572   2.53
# 7 170525_NaHCO3 L + NaCl U_.dxf     7  1.46    1.22  0.656   2.59
# 8 170525_NaHCO3 L + NaCl U_.dxf     8  1.67    1.22  0.860   3.25
# 9 170525_NaHCO3 L + NaCl U_.dxf     9  1.88    1.22  0.461   2.26
# 10 170525_NaHCO3 L + NaCl U_.dxf    10  2.09    1.21  0.522   2.23


# get vendor info
# entire vendor data table
#msdat %>% iso_get_vendor_data_table() %>% knitr::kable()
ven.tab <- iso_get_vendor_data_table(msdat)
str(ven.tab)
ven.tab

rIntensityAll<-ven.tab$`rIntensity All`
rIntensityAll
# double in 'mVs'[16]> # mVs=milliVolt Second (weber)
#[1] 141054.52 141031.51   3848.65 140982.27 140282.23  15604.47  14847.19  14144.31  13479.41  12848.28  12259.87  11713.48
#[13]  11197.01  10698.87  10254.50 140460.79

Rt<-ven.tab$Rt
Rt
# double in 'mVs'[16]>. # mVs=milliVolt Second (weber)
#[1] 141054.52 141031.51   3848.65 140982.27 140282.23  15604.47  14847.19  14144.31  13479.41  12848.28  12259.87  11713.48
#[13]  11197.01  10698.87  10254.50 140460.79


xaxis <- ven.tab$Rt # Rt= time when peak occurs
yaxis <- ven.tab$`rIntensity All`
plot(xaxis,yaxis,type="l")





# get specific parts and add some file information
msdat %>%
  iso_get_vendor_data_table(
    # select peak number, ret. time, overall intensity and all H delta columns
    select = c(Nr., Rt, area = `rIntensity All`, matches("^d \\d+H")),
    # include the Analysis number fron the file info and rename it to 'run'
    include_file_info = c(run = Analysis)
  ) %>%
  knitr::kable()

all_data <- msdat %>% iso_get_all_data()
all_data
# A tibble: 1 x 7
#file_id                   file_type      file_info        raw_data           standards      resistors      vendor_data_table
#<chr>                     <chr>          <list>           <list>             <list>         <list>         <list>
#  1 170525_NaHCO3 L + NaCl U… continuous_fl… <tibble [1 × 15… <tibble [4,297 × … <tibble [6 × … <tibble [3 × … <tibble [16 × 39…


# export to excel
msdat %>% iso_export_to_excel("msdat_file_export") %>% iso_get_file_info(select = c(-file_root)) %>% knitr::kable()
# data sheets available in the exported data file:
readxl::excel_sheets("msdat_file_export.cf.xlsx")
