function r2 = rsquared(rtrue, rhat)
% r2 = rsquared(rtrue, rhat)


r2 = 1-(sum((rtrue(:)-rhat(:)).^2))/(sum((rtrue(:)-mean(rtrue(:))).^2));