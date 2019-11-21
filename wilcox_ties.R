# The second stratum is related with ties when performing Wilcoxon
# All because we are using a non-parametric rank-based test
groupA <- c(0.756701,   0.826087,   0.8136335,  0.9352227,  0.9001585, 
	        0.9056911,  0.832084,   0.9327731,  0.9250814,  0.928702,
	        0.9550562,  0.9355272,  0.9501779,  0.9402093,  0.8571429,
	        0.9108062,  0.85209,    0.8650794,  0.8185328,  0.8324786,  
	        0.8677686,  0.832981,   0.9278132,  0.9421053)
groupB <- c(0.08415041, 0.09426329, 0.07260758, 0.17096129, 0.14607775,
	        0.10448069, 0.07972549, 0.10642058, 0.07338358, 0.11341897, 
	        0.12265566, 0.10894262, 0.13963892, 0.16418979, 0.10338274,
	        0.12913419, 0.13887399)
wilcox.test(as.numeric(groupA), as.numeric(groupB))
plot(density(groupB), xlim=c(0, 1))
lines(density(groupA))

# Let's change the first value to a close value to the second one
groupA[1] <- 0.82608
wilcox.test(as.numeric(groupA), as.numeric(groupB)) # Same p-value
plot(density(groupB), xlim=c(0, 1))
lines(density(groupA))

# Yet if we assign the second value to the first one...
groupA[1] <- groupA[2]
wilcox.test(as.numeric(groupA), as.numeric(groupB)) # Decreased p-value
plot(density(groupB), xlim=c(0, 1))
lines(density(groupA))

# What if we assign a close value to a third value
groupA[3] <- 0.82608
wilcox.test(as.numeric(groupA), as.numeric(groupB)) # Same p-value from before
plot(density(groupB), xlim=c(0, 1))
lines(density(groupA))

# What if a third value equals the second one? Does the p-value decrease?
groupA[3] <- groupA[2]
wilcox.test(as.numeric(groupA), as.numeric(groupB)) # Decreased p-value
plot(density(groupB), xlim=c(0, 1))
lines(density(groupA))

# A fourth value equals the second one
groupA[4] <- groupA[2]
wilcox.test(as.numeric(groupA), as.numeric(groupB)) # Decreased p-value again
plot(density(groupB), xlim=c(0, 1))
lines(density(groupA))

# Even though the distributions are clearly different, the increasing number of
# ties decreases the significance of the difference between groups
