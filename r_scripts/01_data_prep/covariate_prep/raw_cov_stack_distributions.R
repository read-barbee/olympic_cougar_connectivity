################################ Raw cov stack distributions #################################

library(terra)

cov_stack <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_raw_unscaled_30m_11-30-2023.tif")

#remove categorical layers
cov_stack2 <- subset(cov_stack, c(1:16, 20:34))


#first half
plots <- hist(cov_stack2)

#last half
plots2 <- hist(cov_stack2[[17:31]])

#all at once--too large
#plots <- hist(cov_stack2, maxnl = 31)
