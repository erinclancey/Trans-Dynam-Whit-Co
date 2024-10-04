# EpiEstim in Structured Populations

## GENERAL INFORMATION
Repository for the Mathematica notebooks, R code and Whitman County COVID-19 reports for the manuscript:

Erin Clancey $^{1,\ast}$ and Eric T. Lofgren $^1$, Time-varying reproductive number estimation for practical application in structured populations

1. Current Address: Paul G. Allen School for Global Health, Washington State University, Pullman, WA 99164 USA;

  *Corresponding author; e-mail: erin.clancey@wsu.edu;

EpiEstim is a popular statistical framework designed to produce real-time estimates of the time-varying reproductive number, $\mathcal{R}_t$. However, the methods in  EpiEstim have not been tested in small, non-randomly mixing populations to determine if the resulting $\mathcal{R}_t$ estimates are temporally biased. Thus, we evaluate the the temporal performance of EpiEstim estimates when population structure is present. Then, we demonstrate how to recover temporal accuracy using an approximation with EpiEstim. Following a real-world example of a COVID-19 outbreak in a small university town, we generate simulated case report data from a two-population mechanistic model with an explicit generation interval distribution and expression to compute the true $\mathcal{R}_t$. We then compare the time point when $\mathcal{R}_t<1$ to estimates from EpiEstim to quantify the temporal bias. When population structure is present but not accounted for, estimates from EpiEstim prematurely fall below the critical threshold of one. Aggregating incidence data over weeks compared to daily data does not change the affect of population structure. We show it is possible to recover the correct timing by using the lagging subpopulation outbreak to estimate $\mathcal{R}_t$ within EpiEstim. Since population structure can bias $\mathcal{R}_t$ near the critical threshold of one, EpiEstim should be prudently applied to incidence data from structured populations.

This publication was made possible by cooperative agreement CDC-RFA-FT-23-0069 from the CDC’s Center for Forecasting and Outbreak Analytics. Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the Centers for Disease Control and Prevention.  ETL was also funded by NIH R35GM147013.
