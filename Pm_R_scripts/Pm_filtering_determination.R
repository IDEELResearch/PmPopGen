setwd("/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/")
vcftable <- data.table::fread("Pm_HC_raw.table", header = TRUE)

dd <- density(vcftable$QD, na.rm = TRUE)
plot(dd, t = "l")
polygon(dd, col = rgb(1, 0, 0, 0.5))

QD5 <- vcftable |> subset(QD < 5)

dd5 <- density(QD5$QD)
polygon(dd5, col = rgb(1, 0, 0, 0.5))

QD2.5 <- vcftable |> subset(QD < 2.5)

QD_filtered_prop <- nrow(QD2.5)/nrow(vcftable) #0.063 so filtering out 6.3% of variants (30,370)

#going to implement filters to cut down calculation time and get accurate filtering counts

vcftable <- dplyr::anti_join(vcftable, QD2.5)

fs <- density(vcftable$FS)

plot(fs)
fs100 <- subset(vcftable, FS >= 100)
fs100d <- density(fs100$FS)
plot(fs100d)
fs10 <- subset(vcftable, FS > 10)
fs10d <- density(fs10$FS)
plot(fs10d)
fsb10 <- subset(vcftable, FS <= 10)
fsb10d <- density(fsb10$FS)
plot(fsb10d)

FS_filtered_prop <- nrow(fs10)/nrow(vcftable) #0.067 so filtering out 6.7% of remaining variants(30,159)

vcftable <- dplyr::anti_join(vcftable, fs10)

mq <- density(vcftable$MQ)

plot(mq)

mq60 <- subset(vcftable, MQ < 60)

mq60d <- density(mq60$MQ)

plot(mq60d)
polygon(dd5, col = rgb(1, 0, 0, 0.5))

mq40 <- subset(vcftable, MQ < 40)
mq40d <- density(mq40$MQ)
plot(mq40d)
mq50 <- subset(vcftable, MQ < 50)
mq50d <- density(mq50$MQ)
plot(mq50d)
mq55 <- subset(vcftable, MQ < 55)
mq55d <- density(mq55$MQ)
plot(mq55d)
mq59 <- subset(vcftable, MQ < 59)
mq59d <- density(mq59$MQ)
plot(mq59d)

MQ_filtered_prop <- nrow(mq50)/nrow(vcftable) #0.145 so filtering out 14.5% of remaining variants(60,730)

vcftable <- dplyr::anti_join(vcftable,mq50)

mqrs <- density(vcftable$MQRankSum, na.rm = TRUE)
plot(mqrs)

mqrs5 <- subset(vcftable, MQRankSum <= -5)
mqrs5d <- density(mqrs5$MQRankSum)
plot(mqrs5d)
mqrs2.5 <- subset(vcftable, MQRankSum < -2.5)
mqrs2.5d <- density(mqrs2.5$MQRankSum)
plot(mqrs2.5d)

MQRS_filtered_prop <- nrow(mqrs2.5)/nrow(vcftable) #0.012 so filtering out 1.2% of remaining variants(4,344)

vcftable <- dplyr::anti_join(vcftable,mqrs2.5)

rprs <- density(vcftable$ReadPosRankSum, na.rm = TRUE)

plot(rprs)

plot(rprs, xlim=c(-10,-2))

plot(rprs, xlim=c(-5,10))
plot(rprs, xlim=c(-5,6))

rprs5 <- subset(vcftable, ReadPosRankSum < -5)

rprs5d <- density(rprs5$ReadPosRankSum)

plot(rprs5d)

rprs2.5 <- subset(vcftable, ReadPosRankSum < -2.5)

rprs2.5d <- density(rprs2.5$ReadPosRankSum)

plot(rprs2.5d)

RPRS_filtered_prop <- nrow(rprs2.5)/nrow(vcftable) #0.005 so filtering out 0.5% of remaining variants(1,887)

filteredtable <- read.table("hardfiltered.table", header = TRUE)
biallelictable <- read.table("filtered_biallelics.table", header = TRUE)