# x1 <- c(539,
#         573.5,
#         624.5,
#         628.5,
#         609.5
        
        
# )
# y1 <- c(317,
#         320.5,
#         321.5,
#         313.5,
#         288.5
        
        
# )
# x2 <- c(575,
#         602.5,
#         581.5,
#         559.5,
#         519.5,
#         532.5
        
# )
# y2 <- c(311,
#         285.5,
#         299.5,
#         311.5,
#         303.5,
#         314.5
        
# )


# erg <- matrix(nrow=length(x1), ncol = length(x2))
# for (i in 1:length(x1)){
#   for (j in 1:length(x2)){
#     distance <- sqrt((x1[i]-x2[j])^2 + (y1[i]-y2[j])^2)
#     erg[i,j] <- distance
#   }
# }
# erg

#Der ganze Teil habe ich nur für die DIstanzberechnung benutzt


#######################

differences <- read.table("D:/Uni/Sommersemester_25/z_projekt/Filopodia_Base_Quantification/R3/differences_bases.txt",
                          col.names = "Value", colClasses = "numeric")
differences <- differences/2
differences <- differences/25.2 ##Umrechenung 
h <- hist(differences[,1][differences[,1] > 0], freq=FALSE, breaks = 14)
x <- h$breaks
# y <- h$density
y <- h$counts
x <- x[x<66]

plot(x,y)
lamb_hat <- 1/mean(y)
y_hat<-dexp(x,rate=1/lamb_hat)
plot(x,y_hat)

library(EnvStats)
eexp(y) ### Fitting einer Exponentialverteilung


plot(0:55, sort(rexp(0:55,rate=0.05882353 ),decreasing = TRUE))

plot(x,h$density)
hist(differences[,1][differences[,1] > 0], freq=FALSE, breaks = 14)
lines(x, 0.05882353*exp(x*-0.05882353) )

cor(h$density,0.05882353*exp(x*-0.05882353))^2

##############

#Theoretische Berechnung von Birth and Death rate
#Wahrscheinlich nicht wirklich relevant

h2 <- hist(differences[,1], freq=FALSE, breaks=20)

est <- enorm(differences[,1])

est$parameters['sd']

lines(seq(-2,2,0.01), dnorm(seq(-2,2,0.01), mean=0, sd = est$parameters['sd']))


lengths <- read.table("D:/Uni/Sommersemester_25/z_projekt/Filopodia_Base_Quantification/R3/distances.txt",
                      col.names = "Value", colClasses = "numeric")
lengths <- lengths/25.2
hist(lengths[,1], breaks=25)
mean(lengths[,1])

mean(lengths[,1])

mean_l <- mean(lengths[,1])
death_rate <- est$parameters['sd']

### mean_l = birth_rate/death_rate
### -> mean_l * death_rate = birth_rate
birth_rate <- death_rate*(mean_l*1000/7)
death_rate
birth_rate
cor(h2$density, dnorm(h2$breaks[1:17],mean=0, sd=rate))^2

### in actin 
differences_actin <- differences*1000/7
length_actin <- lengths*1000/7
mean_la <- mean(length_actin[,1])
est_act <- enorm(differences_actin[,1])
hist(differences_actin[,1], breaks=20)
rate_la <- est_act$parameters['sd']
mean_la*rate_la
