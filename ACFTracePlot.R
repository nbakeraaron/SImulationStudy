
##### ACF Plotting ######

Suff.Beta.ACF <- BetaAggregate(Sufficient.BetaChains, "Sufficient")
Anc.Beta.ACF <- BetaAggregate(Ancillary.BetaChains, "Ancillary")
ASIS.Beta.ACF <- BetaAggregate(ASIS.BetaChains, "ASIS")
Sand.Beta.ACF <- BetaAggregate(Sandwich.BetaChains, "Sandwich")

ACF <- rbind(Suff.Beta.ACF, Anc.Beta.ACF, ASIS.Beta.ACF, Sand.Beta.ACF)

ggplot(ACF, aes(x = Lag, y = values, colour = Method)) + geom_line() + facet_grid(Beta ~ Chain)


### Trace Plotting ####

trace(Sufficient.BetaChains, 1, Beta.3)


