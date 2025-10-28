
hammingDistance=function(G1, G2){
  allMistakeOne=FALSE
  if(allMistakeOne){
    Gtmp=(G1+G2)%%2
    Gtmp=Gtmp+t(Gtmp)
    nrRevesals=sum(Gtmp==2)/2
    nrIncDel=sum(Gtmp==1)/2
    hammingDis=nrRevesals+nrIncDel
    
  }else{
    hammingDis=sum(abs(G1-G2))
    hammingDis=hammingDis-0.5*sum(G1*t(G1)*(1-G2)*t(1-G2)+G2*t(G2)*(1-G1)*t(1-G1))
  }
  return(hammingDis)
}


Evaluation.DAG=function(estimated.adjace,true.adjace,K){
  
  p=ncol(true.adjace[[1]])
  ## Converting into 01
  estimated.adjace_01=list()
  for (k in 1:K) {
    estimated.adjace_01[[k]] = estimated.adjace[[k]]
    estimated.adjace_01[[k]][which(estimated.adjace[[k]]!=0)]=1
  }
  
  true.adjace_01=list()
  for (k in 1:K) {
    true.adjace_01[[k]] = true.adjace[[k]]
    true.adjace_01[[k]][which(true.adjace[[k]]!=0)]=1
  }
  
  ### Recall
  Recall = 0
  for (k in 1:K) {
    Recall= Recall + ifelse(sum(true.adjace_01[[k]]), sum(true.adjace_01[[k]]*estimated.adjace_01[[k]])/sum(true.adjace_01[[k]]), 0)
  }
  Recall = Recall/k
  
  ### FDR
  FDR = 0
  for (k in 1:K) {
    FDR = FDR + ifelse(sum(estimated.adjace_01[[k]]), 1 - sum(true.adjace_01[[k]]*estimated.adjace_01[[k]])/sum(estimated.adjace_01[[k]]), 0)
  }
  FDR = FDR/K
  
  ### MCC
  MCC = 0
  for (k in 1:K) {
    MCC = MCC + EvaluationMeasures::EvaluationMeasures.MCC(as.vector(true.adjace_01[[k]]),as.vector(estimated.adjace_01[[k]]))
  }
  MCC = MCC/K
  
  ### HM
  HM = 0
  for (k in 1:K) {
    HM = HM + hammingDistance(estimated.adjace_01[[k]],true.adjace_01[[k]])/(p*(p-1))
  }
  HM = HM/K
  
  
  if (is.nan(MCC)){MCC=0}
  Eval_result=c(Recall,FDR,HM,MCC)

  return(list(Eval_result=Eval_result))
}
