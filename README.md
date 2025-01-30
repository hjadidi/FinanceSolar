# FinanceSolar

Default Probability Estimation via Closed-form Formula and Monte Carlo Simulation:<br>

This code is associated with the published paper in the Energy Economics Journal, titled "Risk Mitigation in Project Finance for Utility-Scale Solar PV Projects". The full paper is available at: https://doi.org/10.1016/j.eneco.2025.108221. The code aims to estimate the cumulative first passage default probability for project finance loans related to utility-scale solar PV projects. This probability is crucial for pricing Credit Default Swaps. Detailed information about the methodology can be found in the mentioned manuscript. <br>

The code's output, represented by "df," provides the cumulative default probability for each year of the loan life. It utilizes two approaches: Monte Carlo Simulation and the Closed-Form formula (Equation 6 in the manuscript), with a specified leverage ratio (L). Users have the flexibility to adjust input data, including the energy yield of a solar power plant and project finance loan characteristics. Additionally, it is essential to select an appropriate uncertainty level for the energy yield, which can be estimated by referring to [1], [2]. <br>
<br>
Ref: <br>
[1] Thevenard, Didier, and Sophie Pelland. "Estimating the uncertainty in long-term photovoltaic yield predictions." Solar energy 91 (2013): 432-445. <br>
[2] Jadidi, Hossein, et al. "Bayesian updating of solar resource data for risk mitigation in project finance." Solar Energy 207 (2020): 1390-1403. <br>
