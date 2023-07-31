#Set working directory (Remember to change from "/" to "/")
library("gridExtra")
library("cowplot")

setwd("C:/Users/afua/S3_bucket/Partitioning_ET2/Instance_Partitioning_Africa/Validation_and_stats_best3_crop")
print(getwd())

data <- read.csv("Cropland_area_csv2.csv", header= TRUE)
head(data)

par(mar = c(5, 5, 1, 1)) # c(bottom, left, top, right)
plot <-    
  plot(data$FAO, data$Wapor, pch=19, col ="black",
       xlab= "FAO cultivated area (" ~ km ^2 ~ ")", 
       ylab= "Cropland area in landcover products (" ~ km ^2 ~ ")",
       cex.lab = 1.1, cex.axis=1.1)

abline(WaPor <- lm(data$Wapor ~ data$FAO), col = "black", lwd = 3, lty = 2 )

points(data$FAO, data$Copernicus, col = "#990033", pch = 15)
abline(Copernicus <- lm(data$Copernicus ~ data$FAO), col = "#990033", lwd = 3, lty = 3)

points(data$FAO, data$Glc_share, col = "#0000CC", pch = 17)
abline(Glc_share <- lm(data$Glc_share ~ data$FAO), col = "#0000CC", lwd = 3, lty = 4)

points(data$FAO, data$Modis, col = "#CC6600", pch = 10)
abline(Modis <- lm(data$Modis ~ data$FAO), col = "#CC6600", lwd = 3, lty = 5)

points(data$FAO, data$Sentinel, col = "#33FF99", pch = 20)
abline(Sentinel <- lm(data$Sentinel ~ data$FAO), col = "#33FF99", lwd = 3, lty = 6)

points(data$FAO, data$CCI, col = "#660066", pch = 21)
abline(CCI <- lm(data$CCI ~ data$FAO), col = "#660066", lwd = 3, lty = 1)

points(data$FAO, data$Globecover, col = "#CCCC66", pch = 18)
abline(Globecover <- lm(data$Globecover ~ data$FAO), col = "#CCCC66", lwd = 3, lty = 1)

#Calculate correlations
format(summary(WaPor)$adj.r.squared, digits=4)
format(summary(Copernicus)$adj.r.squared, digits=4)
format(summary(Glc_share)$adj.r.squared, digits=4)
format(summary(Modis)$adj.r.squared, digits=4)
format(summary(Sentinel)$adj.r.squared, digits=4)
format(summary(CCI )$adj.r.squared, digits=4)
format(summary(Globecover)$adj.r.squared, digits=4)


legend( x = "topleft",
        legend = c("GLC-SHARE" ~ R^2 ~ ": 0.85","WaPOR" ~ R^2 ~ ": 0.82", "Copernicus" ~ R^2 ~ ": 0.81", "CCI" ~ R^2 ~ ": 0.76", "MODIS" ~ R^2 ~ ": 0.72", "GlobCover" ~ R^2 ~ ": 0.51", "Sentinel-2" ~ R^2 ~ ": 0.36"),
        col = c("#0000CC","black","#990033","#660066","#CC6600", "#CCCC66", "#33FF99"),
        lwd=2, lty=c(4, 2, 3, 5, 6, 1,1),
         cex = 1.1, horiz = F )
