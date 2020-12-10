
#Fig. S2a
masterpeaks<- data.frame (read.csv ("MasterPeaks.csv", as.is=T)) #MasterPeaks file is Table S1

library (inflection)

quantile (masterpeaks$Peak.Score, seq (0,1,.1))
#0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
#2.06      5.05      6.21      7.30      8.86     10.90     14.30     20.10     32.80     69.80 986657.75 



#range based on quantiles 
0.2*32247
#6449
0.8*32247
#25798
#trimmed to  6449:25797 (removed top and bottom 20% to get a convex curve to determine the inflection point)


check_curve (6449:25797, sort (masterpeaks$Peak.Score[6449:25797]))
#$ctype
#[1] "convex"

trimmed<- masterpeaks$Peak.Score [6449:25798]
trimmed<- sort (trimmed)
bede (1:19350, trimmed, 0) #0 is the index (if data is convex/concave then index=0)
#$iplast
#[1] 16153  

#$iters
#n a     b   EDE
#1 19350 1 19350 16153.

masterpeaks$Peak.Score[16153]
#[1] 10.9  #inflection point corresponded to a peak score of 10.9- set as cutoff.








