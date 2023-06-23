

#DTPP evaluating drug-direct target interaction by DDS which is determined by logit model by combining DETP socre (DES) and DBTP model score (BS) 

model<-glm(label~DES+BS,data = train_data, family = binomial(link ="logit"))
#
summary(model_1)
#prediction
fitted.prob<-predict(model, newdata = test_data, type = "response")
test_data$DDS<-model_1$fitted.values
 


