# NOTE: Application I starts at row 362, Application II starts at row 478
########################## BACHELOR THESIS ON BINOMIAL TREES IN OPTION PRICING ##########################


############# IMPLEMENTING THE FUNCTIONS #############

# static methods
evaluateCallPayoff <- function(S, K) pmax(S - K, 0)
evaluatePutPayoff <- function(S, K) pmax(K - S, 0)

.parseOptionType <- function(optionType) {
  # Parses optionType (string length 1) and 
  # raises possible errors in the input
  # otherwise always returns the correct type 
  # in a length 4 vector in the order: 
  # binary/vanilla, long/short, american/european, call/put
  # when missing values vanilla, long, american are presumed
  # typos or additional words are not checked
  
  input <- strsplit(tolower(optionType), " ")[[1]]
  isLong <- FALSE
  isVanilla <- FALSE
  isEuropean <- FALSE
  isShort <- FALSE
  isBinary <- FALSE
  isAmerican <- FALSE
  isCall <- FALSE
  isPut <- FALSE
  
  if ("short" %in% input) isShort <- TRUE           else isLong <- TRUE
  if ("binary" %in% input) isBinary <- TRUE         else isVanilla <- TRUE
  if ("european" %in% input) isEuropean <- TRUE  
  else if ("american" %in% input)isAmerican <- TRUE else {cat("American option assumed\n"); isAmerican <- TRUE}
  
  if ("call" %in% input)        isCall <- TRUE
  if ("put" %in% input)         isPut <- TRUE
  if ("long" %in% input)        isLong <- TRUE
  if ("vanilla" %in% input)     isVanilla <- TRUE
  if ("american" %in% input)    isAmerican <- TRUE
  if (isLong && isShort)        stop("optionType cannot be both long and short")
  if (isVanilla && isBinary)    stop("optionType cannot be both vanilla and binary")
  if (isEuropean && isAmerican) stop("optionType cannot be both European and American")
  if (isCall && isPut)          stop("optionType cannot be both call and put")
  if (!(isCall || isPut))       stop("call/put not specified")
  
  
  res <- c("long", "short", "vanilla", "binary", "european", "american", "call", 
           "put")[c(isLong, isShort, isVanilla, isBinary, isEuropean, isAmerican, isCall, isPut)]
  return(res)
}


priceCRR <- function(N, S0, U, D, er, K, optionType, calculateExerciseDates=TRUE, dividend=0) {
  N <- N + 1 # Adding timestep 0
  rnProb <- (er - D) / (U - D)
  
  type <- if (length(optionType) == 4) optionType else .parseOptionType(optionType)
  
  numberOfNodes <- (N * (N + 1) / 2)
  parentNodes <- 1:(N * (N - 1) / 2)
  endNodes <- (N * (N - 1) / 2 + 1):numberOfNodes
  connections <- {
    cbind(parentNodes, 
          parentNodes + rep(1:(N-1), 1:(N-1)), 
          parentNodes + rep(1:(N-1), 1:(N-1)) + 1)
  }
  
  ########## PRICING THE TREE #########
  priceTree <- function() {
    
    evaluatePayoff <- {
      if (type[2] == "binary" && type[4] == "call") function(S, K) as.numeric(S >= K)
      else if (type[2] == "binary") function(S, K) as.numeric(S <= K)
      else if (type[4] == "call") evaluateCallPayoff
      else evaluatePutPayoff
    }
    
    valueEuropean <- function(nodeNumber) {
      if (nodeNumber >= endNodes[1]) {
        return(evaluatePayoff(resultVector["S", nodeNumber], K))
      } else {
        childNodes <- connections[nodeNumber, 2:3]
        return(sum(c(rnProb, 1 - rnProb) * resultVector["Fn", childNodes]) / er)
      }
    }
    
    valueAmerican <- function(nodeNumber) {
      if (nodeNumber >= endNodes[1]) {
        evaluatePayoff(resultVector["S", nodeNumber], K)
      } else {
        noExercise <- valueEuropean(nodeNumber)
        exercise <- evaluatePayoff(resultVector["S", nodeNumber], K)
        
        max(exercise, noExercise)
      }
    }
    
    # checks if premature exercise could be optimal (american but not vanilla call on stock with no dividends)
    isPrematureOptimal <- type[3] == "american" && !(type[2] == "vanilla" && type[4] == "call" && dividend == 0)
    valueNode <- if (isPrematureOptimal) valueAmerican else valueEuropean
    
    
    if (calculateExerciseDates) {
      resultVector <- matrix(c(S0, rep(NA, numberOfNodes * 4 - 1)), 
                             ncol = numberOfNodes, nrow = 4)
      rownames(resultVector) <- c("S", "Fn", "Exercise", "timeStep")
    } else {
      resultVector <- matrix(c(S0, rep(NA, numberOfNodes * 3 - 1)), 
                             ncol = numberOfNodes, nrow = 3)
      rownames(resultVector) <- c("S", "Fn", "timeStep")
    }
    
    resultVector["timeStep", ] <- steps <- rep(1:N, 1:N) - 1
    
    # evaluating values of S
    downMoves <- sequence(1:(N)) - 1
    resultVector["S",] <- S0 * (U^(steps - downMoves) * D^downMoves)
    if (dividend != 0) {
      resultVector["S", endNodes] <- resultVector["S", endNodes] - dividend
    }
    
    # evaluating values of Fn
    for (node in numberOfNodes:1) {
      resultVector["Fn", node] <- valueNode(node)
    }
    
    
    if (calculateExerciseDates) {
      if (!isPrematureOptimal) {
        resultVector["Exercise", parentNodes] <- FALSE
        resultVector["Exercise", endNodes] <- resultVector["Fn", endNodes] > 0
      } else {
        Svalues <- resultVector["S",]
        optionValues <- resultVector["Fn",]
        
        resultVector["Exercise", resultVector["Fn",] == 0] <- FALSE
        
        exerciseNotOptimal <- evaluatePayoff(Svalues, K) != optionValues
        # case where early exercise is optimal at first node
        if ((!exerciseNotOptimal)[1] == TRUE) {
          resultVector["Exercise", 1] <- TRUE
          resultVector["Exercise", -1] <- FALSE
        } else if (sum(exerciseNotOptimal) == length(parentNodes) &&
                   which(exerciseNotOptimal) == parentNodes) {
          # if exercise is only optimal at the end nodes
          resultVector["Exercise", endNodes] <- resultVector["Fn", endNodes] > 0
        }
        
        resultVector["Exercise", exerciseNotOptimal] <- FALSE
        
        # all nodes where premature exercise is optimal
        toConsider <- is.na(resultVector["Exercise",])
        while (any(toConsider)) {
          nodeNumber <- which(toConsider)[1]
          
          if (nodeNumber >= endNodes[1]) { # if none under consideration are parent nodes
            resultVector["Exercise", toConsider] <- resultVector["Fn", toConsider] > 0
            break
          }
          node <- resultVector[,nodeNumber]
          
          stepsRemaining <- N - 1 - node["timeStep"] # at least 1
          
          resultVector["Exercise", nodeNumber] <- TRUE
          # setting subsequent moves only accessible from this node's exercise value to false 
          isPut <- type[4] == "put"
          
          for (i in 1:stepsRemaining) {
            nextNode <- connections[nodeNumber, 2 + isPut] # if call upNode, else downNode
            resultVector["Exercise", nextNode] <- FALSE
            nodeNumber <- nextNode
          }
          toConsider <- is.na(resultVector["Exercise",])
        }
      }
    }
    
    if (type[1] == "short") {
      resultVector["Fn", node] <- (-1) * resultVector["Fn", node]
    }
    
    return(resultVector)
  }
  
  ######## MAIN BODY ########
  values <- priceTree()
  
  # Write what the function returns
  result <- list(optionType = paste(type, collapse = " "), RNProb = rnProb, 
                 value = values, fairPrice = values["Fn", 1])
  return(result)
}

drawTreeGraph <- function(resultMatrix, values, title, 
                          exerciseData = FALSE, toScale = FALSE, 
                          showValues = TRUE) {  
  N <- resultMatrix["timeStep", ncol(resultMatrix)] + 1
  
  numberOfNodes <- (N * (N + 1) / 2)
  parentNodes <- 1:(N * (N - 1) / 2)
  endNodes <- (N * (N - 1) / 2 + 1):numberOfNodes
  connections <- {
    cbind(parentNodes, 
          parentNodes + rep(1:(N-1), 1:(N-1)), 
          parentNodes + rep(1:(N-1), 1:(N-1)) + 1)
  }
  
  u <- 1
  d <- -1
  
  downMoves <- sequence(1:N) - 1
  steps <- rep(1:N, 1:N) - 1
  
  # If toScale, values of first given value is used for y position
  yPos <- if (toScale) resultMatrix[values[1],] 
          else u*(steps - downMoves) + d*downMoves
  xPos <- resultMatrix["timeStep",]
  
  endNodeRange <- yPos[endNodes][c(length(endNodes), 1)]
  endSize <- max(endNodeRange) - min(endNodeRange)
  endMiddle <- mean(endNodeRange)
  newSize <- endSize*1.1
  newNodeRange <- endMiddle + c(-newSize / 2, newSize / 2)
  
  plot(1, 1, col="white", main = title, 
       xlab = "Time step", 
       ylab = ifelse(toScale, paste("Value of", values), ""), 
       xlim = c(0, N),
       ylim = newNodeRange)
  
  
  upNodes <- connections[,2]
  downNodes <- connections[,3]
  segments(xPos[parentNodes], yPos[parentNodes], 
           xPos[upNodes], yPos[upNodes], 
           col="cyan3", lwd=2)
  
  segments(xPos[parentNodes], yPos[parentNodes], 
           xPos[downNodes], yPos[downNodes], 
           col="cyan3", lwd=2)
  
  points(xPos, yPos, pch = 19, lwd = 6, col = "cyan3")
  
  
  if (exerciseData) {
    toDraw <- as.logical(resultMatrix["Exercise",])
    points(xPos[toDraw], yPos[toDraw],
           pch = 19, lwd = 6, col = "red")
  }
  
  colorsUsed <- c("darkorange",
                  "magenta3",
                  "turquoise4",
                  "green3",
                  "darkorange3",
                  "purple",
                  "brown",
                  "deeppink",
                  "darkturquoise",
                  "darkolivegreen",
                  "midnightblue"
  )[1:length(values)]
  
  fontSize <- 0.9
  textXPos <- 0.5
  textYPos <- -0.7
  
  if (showValues) {
    for (i in seq_along(values)) {
      text(xPos, yPos, round(resultMatrix[values[i],], 2), 
           adj = c(textXPos, textYPos), col=colorsUsed[i], 
           cex = fontSize)
      textYPos <- textYPos + 2.5
    }
  }
  
  legend("topleft", legend = values, fill = colorsUsed, 
         bg = "transparent", box.lty = 0)
  if (exerciseData) {
    legend("topleft", 
           legend = c(values, "Points of optimal exercise"), 
           fill = c(colorsUsed, "red"), 
           box.lty = 0, bg = "transparent")
  }
}

valueBSM <- function(S0, K, mu, sigma, er, N, optionType) {
  type <- if (length(optionType) == 4) optionType 
          else .parseOptionType(optionType)
  if (type[3] == "american") 
    stop("This function can only value european options")
  
  d1 <- (log(S0/K) + (log(er) + sigma^2 / 2)*N) / (sigma*sqrt(N))
  d2 <- d1 - sigma * sqrt(N)
  
  multiplier <- ifelse(type[1] == "short", -1, 1)
  
  if (type[2] == "binary" && type[4] == "call") 
    (1 / er^N)  * pnorm(d2) * multiplier
  else if (type[2] == "binary") 
    (1 / er^N)  * pnorm(-d2) * multiplier
  else if (type[4] == "call") 
    (pnorm(d1) * S0 - (1 / er^N) * K * pnorm(d2)) * multiplier
  else 
    ((1 / er^N) * K * pnorm(-d2) - pnorm(-d1) * S0 ) * multiplier
}

# Only relevant to Application I
generateTrackingPortfolio <- function(resultsMatrix, u, d, er) {
  N <- resultsMatrix["timeStep", ncol(resultsMatrix)]
  N <- N + 1
  numberOfNodes <- (N * (N + 1) / 2)
  parentNodes <- 1:(N * (N - 1) / 2)
  endNodes <- (N * (N - 1) / 2 + 1):numberOfNodes
  connections <- {
    cbind(parentNodes, 
          parentNodes + rep(1:(N-1), 1:(N-1)), 
          parentNodes + rep(1:(N-1), 1:(N-1)) + 1)
  }
  
  delta <- B <- c(rep(NA, length(parentNodes)), 
                  rep(0, length(endNodes)))
  result <- rbind(resultsMatrix, delta, B)
  
  
  uNodesFn <- resultsMatrix["Fn", connections[, 2]]
  dNodesFn <- resultsMatrix["Fn", connections[, 3]]
  
  result["delta", parentNodes] <- (uNodesFn - dNodesFn) /
    (resultsMatrix["S", parentNodes] * (u - d))
  result["B", parentNodes] <- 
    (u * dNodesFn - d * uNodesFn) / (er * (u - d))
  
  return(result)
}


# Only relevant to Application II
CRRapproximateBSM <- function(N, totalSteps, optionType, prob, S0, K, er, U, D) {
  mu <- prob * (U-1) + (1 - prob) * (D-1)
  sigma <- sqrt(prob * ((U-1) - mu)^2 + (1-prob) * ((D-1) - mu)^2)
  
  Unew <- exp(sigma * sqrt(N / totalSteps))
  Dnew <- Unew^(-1)
  erNew <- er^(N / totalSteps)
  q <- 0.5 + 0.5 * (mu / sigma) * sqrt(N/totalSteps)
  if (q > 1) {
    stop("q greater than 1, q=", q, "\nPlease choose higher number of steps")
  }
  
  mu_hat <- q * log(Unew/Dnew) + log(Dnew)
  sigma_hat <- sqrt(q * (1-q) * log(Unew/Dnew)^2)
  
  fairPrice <- priceCRR(totalSteps, S0, Unew, Dnew, erNew, K, optionType, 
                        calculateExerciseDates = FALSE)$fairPrice
  res <- c(fairPrice, mu_hat, sigma_hat, totalSteps)
  names(res) <- c("fairPrice", "mu_hat", "sigma_hat", "totalSteps")
  
  return(res)
}


################################ Application I ################################


# MS announced a dividend increase to 85 cents per share in their
# Q2 earnings report published at the end of June 2023.
# The ex-dividend date would be on July 28th, which was 
# also the expiration date of some options on the stock
# The dividend would be paid out on august 15th, 
# 12 business days after the ex-dividend date
# The binomial tree will be constructed for the 2 week period
# (10 days) leading up to the ex-dividend date

# rf-rate
# Source: Yahoo Finance
YTBillData <- read.csv("1YTBillData.csv", sep = ",")
YTBillData[,1] <- as.Date(YTBillData[,1])
correctRow <- which(YTBillData[,1] == "2023-07-14")
rfRate <- YTBillData[correctRow, "Close"]
rfRate
rfRate <- as.numeric(substring(rfRate,0,nchar(rfRate)-1))/100 + 1
# This is erT for T = 1
# In the model, one time step = 1 / 252, of a year, 
# since there are on average 252 business days in a year
# hence we have
oneDayRfRate <- rfRate^(1 / 252)
# for 


# Source: optionistics.com
optionPriceHistory <- read.csv("MS options K=91.csv", sep = "\t")
optionPriceHistory
callPriceHistory <- optionPriceHistory[,"callPrice"]


# Stock price data from 13th of July 2018 to 28th of August 2023
stockData <- read.csv("MS.csv", sep = ",")
stockData[,1] <- as.Date(stockData[,1])

head(stockData)
tail(stockData)

plot(stockData$Date, stockData$Adj.Close, type='l', 
     ylab="Adjusted closing prices", xlab="Date")
grid()

fiveYears <- which(as.logical(stockData[,1] < 
                                as.Date("2023-07-14")))
fiveYearPriceHistory <- stockData[fiveYears, "Adj.Close"]

historicNumberOfDays <- length(fiveYearPriceHistory)

historicReturns <- fiveYearPriceHistory[2:historicNumberOfDays] / 
  fiveYearPriceHistory[1:(historicNumberOfDays-1)] - 1
historicMu <- mean(historicReturns)
historicSigma <- sqrt(var(historicReturns))

historicMu
historicSigma


daysToObserve <- which(as.logical(stockData[,1] >= 
                                    as.Date("2023-07-14")))
toObserve <- stockData[daysToObserve,c("Date", "Close")]

nrow(toObserve)

stockValues <- toObserve[,"Close"]
dailyStockReturns <- stockValues[2:11] / stockValues[1:10] - 1

{
  N <- 10 # trading days between 14th and 28th
  S <- S0 <- stockValues[1]
  U <- exp(historicSigma)
  D <- U^(-1)
  er <- oneDayRfRate
  K <- 91
  dividend <- 0.85*er^(-12) # discounted value of dividend
  optionType <- "american call"
}

tree1 <- priceCRR(N, S0, U, D, er, K, optionType, 
                  dividend = dividend)
tpf1 <- generateTrackingPortfolio(tree1$value, U, D, er)
drawTreeGraph(tree1$value, c("S"), 
              "Binomial Tree for Morgan Stanley stock", 
              exerciseData = TRUE, toScale = T)


# On day 1 we write a call option & buy the tracking portfolio

deltas <- Bs <- rep(NA, 10)
fairPrices <- rep(NA, 11)
exerciseData <- rep(NA, 10)


# From day 0 to 9 (1 day before expiration)
for (day in 1:10) {
  tree <- priceCRR(N, S, U, D, er, K, optionType,
                   dividend = dividend)
  fairPrices[day] <- tree$fairPrice
  exerciseData[day] <- tree$value["Exercise", 1]
  
  trackingPf <- 
    generateTrackingPortfolio(tree$value, 
                              U, D, er)[c("delta", "B"),1]
  deltas[day] <- delta <- trackingPf["delta"]
  Bs[day] <- B <- trackingPf["B"]
  
  N <- N - 1
  S <- S * (1 + dailyStockReturns[day])
}

which(as.logical(exerciseData))

# Optimal tracking portfolio value for each day = fair price 
# of option on that day
trackingPortfolios <- deltas * stockValues[1:10] + Bs
trackingPortfolios



# What they were worth the next day (n=1:10)
dayAfterTpf <- deltas * stockValues[2:11] + Bs * er


exerciseDay9 <- evaluateCallPayoff(stockValues[10], K)
exerciseDay10 <- evaluateCallPayoff(stockValues[11], K)

fairPrices[11] <- exerciseDay10

# hedging errors for last days:
hedgingErrorDay9 <- dayAfterTpf[9] - exerciseDay9
hedgingErrorDay10 <- dayAfterTpf[10] - exerciseDay10


PnLOfHedgeExerciseDay9 <- c(dayAfterTpf[1:8] - 
                              trackingPortfolios[2:9],
                            hedgingErrorDay9)
PnLOfHedgeExerciseDay10 <- c(dayAfterTpf[1:9] - 
                               trackingPortfolios[2:10], 
                             hedgingErrorDay10)


# PnL value as of day 0
totalPnLExerciseDay9 <- sum(PnLOfHedgeExerciseDay9*er^(-(1:9)))
totalPnLExerciseDay10 <- sum(PnLOfHedgeExerciseDay10*er^(-(1:10)))

totalHedgingErrorExerciseDay9 <- 
  sum(abs((PnLOfHedgeExerciseDay9)))
totalHedgingErrorExerciseDay10 <- 
  sum(abs((PnLOfHedgeExerciseDay10)))



nearestPath <- c(1, 2, 4, 7, 11, 16, 23, 30, 39, 49, 59)
tree1$value["S",nearestPath]
# Tree from day 1 vs actual stock
drawTreeGraph(tree1$value, c("S"), 
              "Binomial Tree for Morgan Stanley stock", 
              exerciseData = TRUE, toScale = T)
lines(0:10, tree1$value["S", nearestPath], lwd=3, col="darkorchid3")
lines(0:10, stockValues[1:11], lwd=2, col="green3")
legend("topleft", legend = c("S", "Points of optimal exercise", 
                             "Actual stock movement", 
                             "Nearest path"), 
       fill = c("turquoise4", "red", "green3", "darkorchid3"), 
       box.lty = 1, bg = "white")



# Value of option
drawTreeGraph(tree1$value, c("Fn"), 
              "Binomial Tree for Morgan Stanley stock", 
              exerciseData = TRUE, toScale = T, showValues = F)
lines(0:10, fairPrices, lwd=2, col="goldenrod4")
lines(0:10, callPriceHistory, lwd=2, col="forestgreen")
lines(0:10, tree1$value["Fn", nearestPath], lwd=3, col="darkorchid3")
legend("topleft", legend = c("Value of tracking portfolio on each day", 
                  "Historical call prices on each day", 
                  "Points of optimal exercise"), 
       fill = c("goldenrod4", "forestgreen", "red"), 
       box.lty = 1, bg = "white")



summary(lm(callPriceHistory ~ tree1$value["Fn", nearestPath]))

linReg <- lm(callPriceHistory ~ fairPrices)


linReg <- lm(callPriceHistory ~ fairPrices)

qqplot(fairPrices, callPriceHistory, pch=16)
abline(linReg$coefficients)


tree1$fairPrice
stockValues[1:10]

# Hedging error visualised
plot(1:10, PnLOfHedgeExerciseDay10, col="turquoise", pch = 16)
abline(h=0, lty="dashed")
points(9, PnLOfHedgeExerciseDay9[9], col = "red", pch=17)
text(1:10, PnLOfHedgeExerciseDay10, 
     round(PnLOfHedgeExerciseDay10, 4), adj=c(0.5,-0.5), cex = 0.7)
text(9, PnLOfHedgeExerciseDay9[9], 
     round(PnLOfHedgeExerciseDay9[9], 4), adj=c(0.5,-0.5), cex = 0.7)




################################ Application II ################################


{
  N <- 10
  totalSteps <- 1000
  S0 <- 100
  U <- 1.155
  D <- 0.9975
  er <- 1.05
  K <- 200
  prob <- 0.75
  dividend <- 0
  optionType <- "american binary call"
}

myTree <- priceCRR(N, S0, U, D, er, K, optionType, dividend = 0)
riskNeutralProb <- myTree$RNProb

fPamerican <- myTree$fairPrice
drawTreeGraph(myTree$value, "S", 
                "Price movements of underlying", exerciseData = F, toScale = T)

drawTreeGraph(myTree$value, c("Fn", "S"), "American binary call option prices", 
                exerciseData = T, toScale = F)

optionType <- "binary european call"
myTree <- priceCRR(N, S0, U, D, er, K, optionType, dividend = 0)
fPeuropean <- myTree$fairPrice

drawTreeGraph(myTree$value, c("Fn", "S"), "European binary call option prices",
                exerciseData = T, toScale = F)

mu <- riskNeutralProb * (U-1) + (1 - riskNeutralProb) * (D-1)
sigma <- sqrt(riskNeutralProb * ((U-1) - mu)^2 + (1-riskNeutralProb) * ((D-1) - mu)^2)
BSMPrice <- valueBSM(S0, K, mu, sigma, er, N, optionType)

steps <- c(seq(50, 200, 50), seq(250, 2500, 250))

apps <- matrix(rep(NA, length(steps) * 4), nrow = 4)
rownames(apps) <- c("fairPrice", "mu_hat", "sigma_hat", "totalSteps")

# This step takes a while
for (i in seq_along(steps)) {
  apps[,i] <- 
    CRRapproximateBSM(N, steps[i], optionType, riskNeutralProb, S0, K, er, U, D)
  gc(verbose = FALSE, reset = FALSE, full = TRUE)
}


CRR2000 <- priceCRR(2000, 100, U, D, er, K, "european call", F, 0)
qqnorm(log(CRR2000$value["S", CRR2000$value["timeStep",] == 2000]))
qqline(log(CRR2000$value["S", CRR2000$value["timeStep",] == 2000]))


cat(paste("Expected total return", round(mu*N, 7), 
          "\nTotal variance", round(sigma^2*N, 7), 
          "\nFair price", round(BSMPrice, 7), "\n\n"))

appMus <- appVars <- appPrices <- rep(NA, length(steps))

for (i in 1:ncol(apps)) {
  appPrices[i] <- appPrice <- apps["fairPrice", i]
  appMus[i] <- appMu <- apps["mu_hat", i] * steps[i]
  appVars[i] <- appVariance <- apps["sigma_hat", i]^2 * steps[i]
  
  priceDiff <- appPrice - BSMPrice
  muDiff <- appMu - mu*N
  varDiff <- appVariance - sigma^2*N
  
  cat(paste("CRR approximation with", steps[i], "steps", 
            "\nExpected total return", round(appMu, 7), 
            "\tOff by", round(muDiff, 7), "from actual value",
            "\nTotal variance", round(appVariance, 7), 
            "\tOff by", round(varDiff, 7), "from actual value",
            "\nFair price", round(appPrice, 7), 
            "\t\tOff by", round(priceDiff, 7), "from BSM value", "\n\n"))
}

rm(appPrice, muDiff, varDiff)


par(mfrow = c(2, 2))
plot(log10(steps), appMus, pch=16, col="green", main = "Approximations of mu")
abline(h=mu*N, col="darkgreen", lty = "longdash", lwd=2)
points(log10(steps), appMus, pch=16, col="green")

plot(steps, appVars, pch=16, col="orange", main = "Approximations of the variance")
abline(h=sigma^2*N, col="darkorange3", lty = "longdash", lwd = 2)
points(steps, appVars, pch=16, col="orange")

plot(log10(steps), appPrices, pch=16, col="turquoise",
        main = "Approximations of the fair price")
abline(h = BSMPrice, col="blue", lty="longdash", lwd=2)
points(log10(steps), appPrices, pch=16, col="turquoise")

asd <- boxplot((appPrices), horizontal = T, plot = FALSE)
qqnorm(appPrices)
qqline(appPrices)


# time complexity is O(n^2)
system.time(CRRapproximateBSM(N, 1000, optionType, 
                              riskNeutralProb, S0, K, er, U, D))
gc(verbose = FALSE, reset = FALSE, full = TRUE)
system.time(CRRapproximateBSM(N, 2000, optionType, 
                              riskNeutralProb, S0, K, er, U, D))
gc(verbose = FALSE, reset = FALSE, full = TRUE)
system.time(CRRapproximateBSM(N, 4000, optionType, 
                              riskNeutralProb, S0, K, er, U, D))
gc(verbose = FALSE, reset = FALSE, full = TRUE)
# system.time(CRRapproximateBSM(N, 8000, optionType, 
#                               riskNeutralProb, S0, K, er, U, D))

par(mfrow = c(1, 1))

