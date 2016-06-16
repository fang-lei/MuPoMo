
![http://quantnet.wiwi.hu-berlin.de/style/banner.png](http://quantnet.wiwi.hu-berlin.de/style/banner.png)

## ![qlogo](http://quantnet.wiwi.hu-berlin.de/graphics/quantlogo.png) **MuPoMo_data**


```yaml

Name of QuantLet :  MuPoMo_data

Published in :      MuPoMo (Mortality Model for Multip-populations: A Semiparametric Comparison Approach)

Description :      'Reads data from Human Mortality Database and calculate kt, and will be called
                    in MuPoMo_main_twopop and MuPoMo_main_multipop.’

Keywords :         ‘time series, demography, mortality, Lee-Carter method, population, Human Mortality Database’

See also :         ‘MuPoMo_optimization, MuPoMo_normalization, MuPoMo_referencecurve, 
                    MuPoMo_main_twopop, MuPoMo_main_multipop, MuPoMo_bootstrap’

Author :            Lei Fang

Submitted :         Mon, June 13 2016 by Lei Fang


```R


# read multi-pop female mortality rates from Human Mortality Database
shortnames.all = c("AUS", "AUT", "BLR", "BGR", "CAN", "CHL", "CZE", "DNK", "EST", "FIN", "FRATNP", "DEUTNP", "HUN", "ISL", "IRL", 
    "ISR", "ITA", "JPN", "LVA", "LTU", "LUX", "NLD", "NZL_NP", "NOR", "POL", "PRT", "RUS", "SVK", "SVN", "ESP", "CHE", "TWN", "GBR_NP", 
    "USA", "SWE")
names.all = c("Australia", "Austria", "Belarus", "Bulgaria", "Canada", "Chile", "CzechRepublic", "Denmark", "Estonia", "Finland", 
    "France", "Germany", "Hungary", "Iceland", "Ireland", "Israel", "Italy", "Japan", "Latvia", "Lithuania", "Luxembourg", "Netherlands", 
    "NewZealand", "Norway", "Poland", "Portugal", "Russia", "Slovakia", "Slovenia", "Spain", "Switzerland", "Taiwan", "UnitedKingdom", 
    "USA", "Sweden")

loop.all = length(names.all)
for (i in 1:loop.all) {
    nam1 = paste(names.all[i])
    # read data from Human Mortality Database (HMD: http://www.mortality.org/)
    # you need register on HMD to access a free account accound and set up your password
    assign(nam1, hmd.mx(shortnames.all[i], "account", "password", names.all[i]))
    temp1 = hmd.mx(shortnames.all[i], "account", "password", names.all[i])
    nam2  = paste(names.all[i], "lca.female", sep = ".")
    # Lee-Carter method
    assign(nam2, lca(temp1, series = "female", adjust = "dt", interpolate = TRUE))
    temp2 = lca(temp1, series = "female", adjust = "dt", interpolate = TRUE)
    nam3  = paste("ax", names.all[i], "female", sep = ".")  # ax
    assign(nam3, temp2$ax)
    nam4  = paste("bx", names.all[i], "female", sep = ".")  # bx
    assign(nam4, temp2$bx)
    nam5  = paste("kt", names.all[i], "female", sep = ".")  # kt
    assign(nam5, temp2$kt)
}
```
