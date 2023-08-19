### Contents of `wheat.csv`

- `Year`: serves as primary key
- `Q_production`: global quantity of wheat produced, in millions of metric tonnes. Sourced from the FAO data pool for [Crops and livestock products](https://www.fao.org/faostat/en/#data/QCL) in August 2022.
- `Q_demand`: total quantity of wheat consumed, in millions of metric tonnes. Calculated using production from the FAO data pool (see above) and change in wheat stocks sourced from the [USDA Wheat Data page](https://www.ers.usda.gov/webdocs/DataFiles/54282/Wheat%20Data-All%20Years.xlsx?v=8188.3) in August 2022.
- `Q_Ukraine`: quantity of wheat produced in Ukraine, in millions of metric tonnes. Sourced from the FAO data pool for [Crops and livestock products](https://www.fao.org/faostat/en/#data/QCL) in August 2022.
- `P`: global wheat prices, in real 2020 USD per metric tonnes. Sourced from the [FRED What price portal](https://fred.stlouisfed.org/series/PWHEAMTUSDA) in August 2022 and deflated using the CPI.
