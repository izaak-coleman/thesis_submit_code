setwd("Documents/Computing/MscComputingsci/thesis/src/toySMuFin/") # get to data


merge_dat = read.csv("merge_sort_SA_construction_benchmarking.csv")
rad_merge_dat = read.csv("radix_merge_SA_construction_benchmarking_data.csv")
para_rad_dat  = read.csv("parallel_radixSA_construction_benchmarking_data.csv")

# simple way to build discrete plots
rad_merge_levels <- factor(rad_merge_dat$data_size) # factor pulls of the catgor. levels
merge_dat_levels <- factor(merge_dat$data_size)
para_rad_levels  <- factor(para_rad_dat$data_size)

par(mfrow=c(1,2))
plot(merge_dat_levels, log2(merge_dat$ms))  # naive merge sort
plot(rad_merge_levels, log2(rad_merge_dat$ms)) # first radix impl
plot(para_rad_levels, log2(para_rad_dat$ms)) # fastest imp

#continuous plots
par(mfrow=c(1,2))
merge_avs = aggregate(ms~data_size,merge_dat,mean)
rad_merge_avs = aggregate(ms~data_size,rad_merge_dat, mean)
para_rad_avs = aggregate(ms~data_size,para_rad_dat,mean)


# milliseconds to seconds
merge_avs[2] <- merge_avs[2]/1000
rad_merge_avs[2] <- rad_merge_avs[2] / 1000
para_rad_avs[2] <- para_rad_avs[2] / 1000


# reads to kiloreads
merge_avs[1] <- merge_avs[1] / 1000
rad_merge_avs[1] <- rad_merge_avs[1] / 1000
para_rad_avs[1] <- para_rad_avs[1] / 1000


plot(merge_avs$data_size, merge_avs$ms, type="line", col="black", ylab="Run time (s)", xlab="Input data (thousands)", main="SA construction algorithm runtimes")
lines(rad_merge_avs$data_size, rad_merge_avs$ms, type="line", col="red")
lines(para_rad_avs$data_size, para_rad_avs$ms, type="line", col="green")
legend(locator(1),c("merge sort","radix merge", "para radix"),pch=c(1,16), col=c("black","red", "green"))

#ratio plots
merge_rad_merge_ratio <- rad_merge_avs[1]  ## get the explanatory
merge_rad_merge_ratio[2] <-merge_avs[2]/rad_merge_avs[2]  # compute the merge/x ratio

merge_para_rad_ratio <- para_rad_avs[1]
merge_para_rad_ratio[2] <- merge_avs[2]/para_rad_avs[2]

merge_merge_ratio <-merge_avs[1]
merge_merge_ratio[2] <-merge_avs[2]/merge_avs[2]
#plot data
plot(merge_para_rad_ratio$data_size, merge_para_rad_ratio$ms, ylim=c(0,20), type="l", col="green", xlab="Input data (thousands)", ylab="Times speed-up vs merge sort", main="Para radix scales to large data")
lines(merge_para_rad_ratio$data_size, merge_merge_ratio$ms, type="l", col="black")
lines(merge_para_rad_ratio$data_size, merge_rad_merge_ratio$ms,type="l", col="red")


# scripts for generating the cost of merge sort

summed_levels_32 <- read.csv("sorted_accum_levels_32.csv")
summed_levels_64 <- read.csv("sorted_accum_levels_64.csv")
summed_levels_128 <- read.csv("sorted_accum_levels_128.csv")
summed_ints_stack_128 <- read.csv("sorted_accum_ints_stack_128.csv")
summed_ints_heap_128 <- read.csv("sorted_accum_ints_heap_128.csv")




plot(summed_levels_128$level, summed_levels_128$ns/1000000000, xlab="Merge sort recursion level", ylab="Time spent merging at each recursion level (s)", type="line", main="The last merge in lexicographic merge sort is the most expensive", col="red", ylim=c(0,2.55))
lines(summed_levels_64$level, summed_levels_64$ns/1000000000, type="line", col="blue")
lines(summed_levels_32$level, summed_levels_32$ns/1000000000, type="line", col="green")
lines(summed_ints_stack_128$level, summed_ints_stack_128$ns/1000000000, type="line", col="yellow")
lines(summed_ints_heap_128$level, summed_ints_heap_128$ns/1000000000, type="line", col="black")
legend(locator(1),c("128 DNA","64 DNA", "32 DNA", "128 int stack", "128 int heap"),pch=c(1,16), col=c("red","blue", "green", "yellow", "black"))
summed_ints_stack_128
summed_ints_heap_128



sorted_accum_levels <- read.csv("sorted_accum_levels_64.csv")
