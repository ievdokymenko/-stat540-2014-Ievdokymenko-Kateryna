gDat <- read.delim("gapminderDataFiveYear.txt")
str(gDat)
library(plyr)
max(subset(gDat, continent=="Africa")$lifeExp)
(maxLeByCont <- ddply(gDat, ~ continent, summarize, maxLifeExp = max(lifeExp)))
(minGdpByCont <- ddply(gDat, ~ continent, summarize, minGdp = min(gdpPercap)))

ddply(gDat, ~continent, summarize, nUniqCountries = length(unique(country)))

ddply(gDat, ~ continent,
function(x) return(c(nUniqCountries = length(unique(x$country)))))

ddply(gDat, ~ continent, summarize,
      minLifeExp = min(lifeExp), maxLifeExp = max(lifeExp),
      medGdpPercap = median(gdpPercap))

jCountry <- "France"  # pick, but do not hard wire, an example
(jDat <- subset(gDat, country == jCountry))  # temporary measure!
library(lattice)
xyplot(lifeExp ~ year, jDat, type = c("p", "r"))  # always plot the data

jFit <- lm(lifeExp ~ year, jDat)
summary(jFit)

(yearMin <- min(gDat$year))
jFit <- lm(lifeExp ~ I(year - yearMin), jDat)
summary(jFit)
class(jFit)
mode(jFit)

jFun <- function(x) {
  estCoefs <- coef(lm(lifeExp ~ I(year - yearMin), x))
  names(estCoefs) <- c("intercept", "slope")
  return(estCoefs)
}
jFun(jDat) 

jCoefs <- ddply(gDat, ~country, jFun)
str(jCoefs)
tail(jCoefs)

library(xtable)

set.seed(916)
foo <- jCoefs[sample(nrow(jCoefs), size = 15), ]
foo <- xtable(foo)
print(foo, type = "html", include.rownames = FALSE)

jCoefs <- ddply(gDat, ~country + continent, jFun)
set.seed(916)
foo <- jCoefs[sample(nrow(jCoefs), size = 15), ]
foo <- arrange(foo, intercept)
## foo <- foo[order(foo$intercept), ] # an uglier non-plyr way
foo <- xtable(foo)
print(foo, type = "html", include.rownames = FALSE)