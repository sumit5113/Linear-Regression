data<-read.table(file.choose(),header=TRUE,sep=',')

data

attach(data)

clmn_scs_nt_scs <-cbind(pass,not_pass)

clmn_scs_nt_scs

logit <- glm(clmn_scs_nt_scs~s_p_h, family = binomial)

summary(logit)

#Residual deviance:  4.8152  on 6  degrees of freedom  -> chi square

qchisq(p = .95,lower.tail = TRUE,df=6) # 4.8152 < 12.59159 , hence significant the residual model is

plot(s_p_h, fitted.values(logit))
points(s_p_h,pass/not_pass,col='red')
warnings()



