

#DTPP evaluating drug-direct target interaction by DDS which is determined by logit model by combining DETP socre (DES) and DBTP model score (DBS) 

model<-glm(label~DES+DBS,data = train_data, family = binomial(link ="logit"))
#
summary(model)
#prediction
fitted.prob<-predict(model, newdata = test_data, type = "response")
test_data$DDS<-model$fitted.values
 


