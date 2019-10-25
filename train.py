import random
import math
import numpy as np
import bisect
from sklearn import svm

line1=input().split(' ')

increment=float(line1[0])
logScale=int(line1[1])
cutoff=float(line1[2])
myType=int(line1[3])
nType=int(line1[4])

random.seed(1)




class Trainer:
    def __init__(self, rs):
        self.rs=rs
    def nPossibleInsertLocations(self):
        return (len(self.rs)-1)*nType
    
    def loadData(self, name1, name2):
        file1=open(name1, "r")
        file2=open(name2, "r")
        line1=file1.readlines()
        line2=file2.readlines()
        file1.close()
        file2.close()
        random.shuffle(line1)
        random.shuffle(line2)
        dists=[]
        self.label=[]
        for d1, d2 in zip(line1, line2):
            dists.append(list(d1.split(' '))[:])
            self.label.append(1)
            dists.append(list(d2.split(' '))[:])
            self.label.append(0)
            
        del line1
        del line2
        
        data=np.zeros([len(dists), len(self.rs)*nType+1])    
        #calculate structure functions
        for j in range(0, len(dists)):
            for jj in range(1, len(dists[j])-1, 2):
                k=float(dists[j][jj+1])
                name=ord(dists[j][jj])-ord('A')
                index=bisect.bisect(self.rs, k)
                if index!=len(self.rs) and index!=0:
                    index2=index+len(self.rs)*name
                    remainder=(k-self.rs[index-1])/(self.rs[index]-self.rs[index-1])
                    data[j, index2]+=(1.0-remainder)
                    data[j, index2+1]+=remainder        
            data[j, 0]=dists[j][0]
        
        if len(data)==0:
            self.mean=np.zeros(0)
            self.stddev=np.zeros(0)
            self.enabled=np.zeros(0)
            self.processedData=np.zeros(0)
            return
        #normalize data
        self.mean=np.zeros(len(data[0]))    
        self.stddev=np.zeros(len(data[0]))    
        self.enabled=np.zeros(len(data[0]))    
        self.processedData=[[] for _ in range(len(data)) ] 
        for i in range(0, len(data[0])):
            sum=0
            sum2=0
            for j in data:
                sum+=j[i]
                sum2+=j[i]*j[i]
            self.mean[i]=sum/len(data)
            self.stddev[i]=math.sqrt(sum2/len(data)-self.mean[i]*self.mean[i])
            if self.mean[i]>1e-9:
                for j in range(0, len(data)):
                    self.processedData[j].append((data[j][i]-self.mean[i])/self.stddev[i])
                self.enabled[i]=1.0    
        print("data processed, size=", len(self.processedData), end=" ")
        
    def reLoadData(self, name1, name2):
        file1=open(name1, "r")
        file2=open(name2, "r")
        line1=file1.readlines()
        line2=file2.readlines()
        file1.close()
        file2.close()
        random.shuffle(line1)
        random.shuffle(line2)
        dists=[]
        self.label=[]
        for d1, d2 in zip(line1, line2):
            dists.append(list(d1.split(' '))[:])
            self.label.append(1)
            dists.append(list(d2.split(' '))[:])
            self.label.append(0)
            
        del line1
        del line2
        
        data=np.zeros([len(dists), len(self.rs)*nType+1])    
        #calculate structure functions
        for j in range(0, len(dists)):
            for jj in range(1, len(dists[j])-1, 2):
                k=float(dists[j][jj+1])
                name=ord(dists[j][jj])-ord('A')
                index=bisect.bisect(self.rs, k)
                if index!=len(self.rs) and index!=0:
                    index2=index+len(self.rs)*name
                    remainder=(k-self.rs[index-1])/(self.rs[index]-self.rs[index-1])
                    data[j, index2]+=(1.0-remainder)
                    data[j, index2+1]+=remainder        
            data[j, 0]=dists[j][0]
        
        #normalize data
        self.processedData=[[] for _ in range(len(data)) ] 
        for i in range(0, len(data[0])):
            if self.enabled[i]:
                for j in range(0, len(data)):
                    self.processedData[j].append((data[j][i]-self.mean[i])/self.stddev[i])
        print("data processed, size=", len(self.processedData))
        
    def Train(self):
        if len(self.processedData)<= 1000:
            print("too few data, cannot train")
            self.clf=svm.SVC(kernel='linear', C=1);
            self.validateAccuracy=0
            return
        TrainData=self.processedData
        TrainLabel=self.label
               
        self.clf=svm.SVC(kernel='linear', C=1);
        self.clf.fit(TrainData, TrainLabel)
        
        self.trainAccuracy = self.clf.score(TrainData, TrainLabel)
        
    def printResult(self):
        print(self.trainAccuracy, self.validateAccuracy)


    def reValidate(self):
        #redo validation using ALL data
        ValidationData=self.processedData
        ValidationLabel=self.label
        self.validateAccuracy = self.clf.score(ValidationData, ValidationLabel)
        return self.validateAccuracy  
            
    def output(self):
         #transform and output the hyperplane
        cumsumEnabled=np.cumsum(self.enabled)
        cept=self.clf.intercept_
        ci=np.zeros(len(self.mean))
        for i in range(len(self.enabled)):
            if self.enabled[i]:
                ci[i]=self.clf.coef_[0,int(cumsumEnabled[i]-1)]/self.stddev[i]
                cept=cept-self.clf.coef_[0,int(cumsumEnabled[i]-1)]/self.stddev[i]*self.mean[i]
    
        ofilename="newHyperplane.txt"
        ofile=open(ofilename, 'w')
        ofile.write(str(cept[0]))
        ofile.write(' ')
        ofile.write(str(ci[0]))
        ofile.write("\n")
        for j in range(0, nType):
            for i in range(len(self.rs)):
                ofile.write(chr(ord('A')+j))
                ofile.write(' ')
                ofile.write(str(ci[i+j*len(self.rs)+1]))
                ofile.write(' ')
                ofile.write(str(self.rs[i]))
                ofile.write("\n")
        ofile.close()

rs=[]
temp=1.0
while temp<cutoff:
	rs.append(temp)
	if logScale:
		temp=temp*(1.0+increment)
	else:
		temp=temp+increment

t=Trainer(rs)

name1="../trainsoft"+str(myType)+".txt"
name2="../trainhard"+str(myType)+".txt"
t.loadData(name1, name2)
t.Train()

name1="../testsoft"+str(myType)+".txt"
name2="../testhard"+str(myType)+".txt"
t.reLoadData(name1, name2)
acc=t.reValidate()

if acc>0:
    f2=open("validateAccuracy.txt", "w")
    f2.write(str(acc))
    f2.write("\n")
    t.output()
    t.printResult()
