function pc=percentcorrectSDT(x,y)
    pc=normcdf(dprime(x,y)/2);