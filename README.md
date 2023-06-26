# Direct-target-prediction
A simple model used for predicting drug-direct target interactions

This work is based on the fact that the effect connections between drugs and their targets 
can be captured by comparing their gene expression profiles.

We can use the gene expression and chemical data to evaluate the effect and binding connections between drugs and targets, respectively,
and then predit drug direct (effect) targets that bind to the drug and are responsible for the biological effects of the drug.

We caculate the drug-effect target interaction score (DES) between drugs and targets 
using gene expression data by a GSEA approach, and caculate the drug-binding target interaction score (DBS) 
between drugs and targets using classical binary models like DeepPurpose.

Then a binary logistic regression is used to combine DES and DBS and caculates the drug-direct target interaction score (DDS)
between drugs and targets. 
