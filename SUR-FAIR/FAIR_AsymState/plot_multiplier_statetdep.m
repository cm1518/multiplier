

msum=squeeze(median(Msum,1))';
msumlow=squeeze(prctile(Msum,lower,1))';
msumhigh=squeeze(prctile(Msum,upper,1))';
msumlow2=squeeze(prctile(Msum,lower2,1))';
msumhigh2=squeeze(prctile(Msum,upper2,1))';

