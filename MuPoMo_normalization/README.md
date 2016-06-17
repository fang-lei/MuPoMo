![http://quantnet.wiwi.hu-berlin.de/style/banner.png](http://quantnet.wiwi.hu-berlin.de/style/banner.png)

## ![qlogo](http://quantnet.wiwi.hu-berlin.de/graphics/quantlogo.png) **MuPoMo_normalization**


```yaml

Name of QuantLet :  MuPoMo_normalization

Published in :      MuPoMo (Mortality Model for Multip-populations: A Semiparametric Comparison Approach)

Description :      'Normalizes optimal shape variation parameters theta estimated from MuPoMo_optimization, and will be 
                    called in MuPoMo_main_multipop.'

Keywords :         'time series, demography, mortality, population, normalization'

See also :         'MuPoMo_data, MuPoMo_optimization, MuPoMo_referencecurve, 
                    MuPoMo_main_twopop, MuPoMo_main_multipop, MuPoMo_bootstrap'

Author :            Lei Fang

Submitted :         Mon, June 13 2016 by Lei Fang


```

```R


normalization  = function(theta) {
    theta.temp = colSums(theta)
    for (i in 1:loop.31) {
        nam10  = paste("normal.theta", names.31[i], 1, sep = ".")
        assign(nam10, loop.31 * theta[i, 1]/theta.temp[1])
        nam11  = paste("normal.theta", names.31[i], 3, sep = ".")
        assign(nam11, loop.31 * theta[i, 3]/theta.temp[3])
        nam12  = paste("normal.theta", names.31[i], 2, sep = ".")
        assign(nam12, theta[i, 2] - theta.temp[2]/loop.31)
        nam14  = paste("normal.theta", names.31[i], sep = ".")
        assign(nam14, c(eval(parse(text = paste("normal.theta", names.31[i], 1, sep = "."))), eval(parse(text = paste("normal.theta", 
            names.31[i], 2, sep = "."))), eval(parse(text = paste("normal.theta", names.31[i], 3, sep = ".")))))
    }
    ll = list()
    for (i in 1:loop.31) {
        ll[[i]]  = c(eval(parse(text = paste("normal.theta", names.31[i], sep = "."))))
    }
    normal.theta = do.call(rbind, ll)
    return(normal.theta)
}

```
