import os
import sys
import subprocess
from multiprocessing import Process

import random
import math
import numpy as np
import bisect
from sklearn import svm




OMP_NUM_THREADS = "20"


INCR = 25000
MAX_DATASET = 25000

D2_MIN_BOUND = 0.05
D2_MIN_RANGE = 0.013
N_RANGE_COEFF = 10.0
PREFIX_A = "train"
PREFIX_B = "test"

FORECAST_SPAN = 50


MINTRAJ_A = 125000 + (FORECAST_SPAN * INCR)
MAXTRAJ_A = 250000000
MINTRAJ_B = MAXTRAJ_A
MAXTRAJ_B = 380150000

PROG = "Debug"

DIRLIST = [1, 2]


increment = 0.1
logScale = 1
cutoff = 10.0


class Trainer:
    def __init__(self, rs, nType):
        self.rs = rs
        self.nType = nType

    def nPossibleInsertLocations(self):
        return (len(self.rs)-1)*self.nType

    def loadData(self, name1, name2):
        file1 = open(name1, "r")
        file2 = open(name2, "r")
        line1 = file1.readlines()
        line2 = file2.readlines()
        file1.close()
        file2.close()
#        random.shuffle(line1)
#        random.shuffle(line2)
        dists = []
        self.ts_fit = []
        self.label = []
        for d1, d2 in zip(line1, line2):
            line_soft = d1.split(' ')
            line_hard = d2.split(' ')
            dists.append(line_soft[1:])
            self.ts_fit.append(line_soft[0])
            self.label.append(1)
            dists.append(line_hard[1:])
            self.ts_fit.append(line_hard[0])
            self.label.append(0)
        with open ("ts_fit.txt", "w") as f:
            for item in self.ts_fit:
                f.write("{}\n".format(item))

        del line1
        del line2

        data = np.zeros([len(dists), len(self.rs)*self.nType+1])
        #calculate structure functions
        for j in range(0, len(dists)):
            for jj in range(1, len(dists[j])-1, 2):
                k = float(dists[j][jj+1])
                name = ord(dists[j][jj])-ord('A')
                index = bisect.bisect(self.rs, k)
                if index != len(self.rs) and index != 0:
                    index2 = index+len(self.rs)*name
                    remainder = (k-self.rs[index-1]) / \
                        (self.rs[index]-self.rs[index-1])
                    data[j, index2] += (1.0-remainder)
                    data[j, index2+1] += remainder
            data[j, 0] = dists[j][0]

        if len(data) == 0:
            self.mean = np.zeros(0)
            self.stddev = np.zeros(0)
            self.enabled = np.zeros(0)
            self.processedData = np.zeros(0)
            return
        #normalize data
        self.mean = np.zeros(len(data[0]))
        self.stddev = np.zeros(len(data[0]))
        self.enabled = np.zeros(len(data[0]))
        self.processedData = [[] for _ in range(len(data))]
        for i in range(0, len(data[0])):
            sum = 0
            sum2 = 0
            for j in data:
                sum += j[i]
                sum2 += j[i]*j[i]
            self.mean[i] = sum/len(data)
            self.stddev[i] = math.sqrt(
                sum2/len(data)-self.mean[i]*self.mean[i])
            if self.mean[i] > 1e-9:
                for j in range(0, len(data)):
                    self.processedData[j].append(
                        (data[j][i]-self.mean[i])/self.stddev[i])
                self.enabled[i] = 1.0
        print("data processed, size=", len(self.processedData))

    def reLoadData(self, name1, name2):
        file1 = open(name1, "r")
        file2 = open(name2, "r")
        line1 = file1.readlines()
        line2 = file2.readlines()
        file1.close()
        file2.close()
        #random.shuffle(line1)
        #random.shuffle(line2)
        dists = []
        self.ts_val = []
        self.label = []
        for d1, d2 in zip(line1, line2):
            line_soft = d1.split(' ')
            line_hard = d2.split(' ')
            dists.append(line_soft[1:])
            self.ts_val.append(line_soft[0])
            self.label.append(1)
            dists.append(line_hard[1:])
            self.ts_val.append(line_hard[0])
            self.label.append(0)
        with open ("ts_val.txt", "w") as f:
            for elem in self.ts_val:
                f.write("{}\n".format(elem))
            


        del line1
        del line2

        data = np.zeros([len(dists), len(self.rs)*self.nType+1])
        #calculate structure functions
        for j in range(0, len(dists)):
            for jj in range(1, len(dists[j])-1, 2):
                k = float(dists[j][jj+1])
                name = ord(dists[j][jj])-ord('A')
                index = bisect.bisect(self.rs, k)
                if index != len(self.rs) and index != 0:
                    index2 = index+len(self.rs)*name
                    remainder = (k-self.rs[index-1]) / \
                        (self.rs[index]-self.rs[index-1])
                    data[j, index2] += (1.0-remainder)
                    data[j, index2+1] += remainder
            data[j, 0] = dists[j][0]

        #normalize data
        self.processedData = [[] for _ in range(len(data))]
        for i in range(0, len(data[0])):
            if self.enabled[i]:
                for j in range(0, len(data)):
                    self.processedData[j].append(
                        (data[j][i]-self.mean[i])/self.stddev[i])
        print("data processed, size=", len(self.processedData))

    def Train(self):
        if len(self.processedData) <= 1000:
            print("too few data, cannot train")
            self.clf = svm.SVC(kernel='linear', C=1)
            self.validateAccuracy = 0
            return
        TrainData = self.processedData
        TrainLabel = self.label

        self.clf = svm.SVC(kernel='linear', C=1)
        self.clf.fit(TrainData, TrainLabel)

        self.trainAccuracy = self.clf.score(TrainData, TrainLabel)

    def printResult(self):
        print(self.trainAccuracy, self.validateAccuracy)

    def reValidate(self):
        #redo validation using ALL data
        self.ValidationData = self.processedData
        self.ValidationLabel = self.label
        self.validateAccuracy = self.clf.score(self.ValidationData, self.ValidationLabel)
        return self.validateAccuracy
    
    def predict(self):
        self.predicted_labels = self.clf.predict(self.ValidationData)
        return self.predicted_labels



    def output(self):
         #transform and output the hyperplane
        cumsumEnabled = np.cumsum(self.enabled)
        cept = self.clf.intercept_
        ci = np.zeros(len(self.mean))
        for i in range(len(self.enabled)):
            if self.enabled[i]:
                ci[i] = self.clf.coef_[
                    0, int(cumsumEnabled[i]-1)]/self.stddev[i]
                cept = cept - \
                    self.clf.coef_[0, int(cumsumEnabled[i]-1)] / \
                    self.stddev[i]*self.mean[i]

        ofilename = "newHyperplane.txt"
        ofile = open(ofilename, 'w')
        ofile.write(str(cept[0]))
        ofile.write(' ')
        ofile.write(str(ci[0]))
        ofile.write("\n")
        for j in range(0, self.nType):
            for i in range(len(self.rs)):
                ofile.write(chr(ord('A')+j))
                ofile.write(' ')
                ofile.write(str(ci[i+j*len(self.rs)+1]))
                ofile.write(' ')
                ofile.write(str(self.rs[i]))
                ofile.write("\n")
        ofile.close()
        


def mytrainer(myType, nType):
    rs = []
    temp = 1.0
    while temp < cutoff:
        rs.append(temp)
        if logScale:
            temp = temp*(1.0+increment)
        else:
            temp = temp+increment

    t = Trainer(rs, nType)
    name1 = "../trainsoft"+str(myType)+".txt"
    name2 = "../trainhard"+str(myType)+".txt"
    t.loadData(name1, name2)
    t.Train()

    name1 = "../testsoft"+str(myType)+".txt"
    name2 = "../testhard"+str(myType)+".txt"
    t.reLoadData(name1, name2)
    acc = t.reValidate()

    if acc > 0:
        f2 = open("validateAccuracy.txt", "w")
        f2.write(str(acc))
        f2.write("\n")
        t.output()
        t.printResult()
    predicted_labels =t.predict()
    with open("ts_validation.txt", "w") as f:
        f.write("ts,predicted_label,truth\n")
        for i,label in enumerate(predicted_labels):
            f.write("{},{},{}\n".format(t.ts_val[i], label, t.ValidationLabel[i]))
        

    


if __name__ == "__main__":
    random.seed(1)
    prog = "../a.out"
    os.putenv("OMP_NUM_THREADS", OMP_NUM_THREADS)
    cwd = os.getcwd()
    for i in DIRLIST:
        path = cwd + "/" + str(i)
        try:
            os.mkdir(path)
        except:
            pass
        os.chdir(path)
        input1 = "{} {} {} {} {} {} {} {} {} {} {} Exit\n".format(PROG, MINTRAJ_A, MAXTRAJ_A,
                                                                  INCR, MAX_DATASET, N_RANGE_COEFF, PREFIX_A, i, D2_MIN_BOUND, D2_MIN_RANGE, FORECAST_SPAN)
        input2 = "{} {} {} {} {} {} {} {} {} {} {} Exit\n".format(PROG, MINTRAJ_B, MAXTRAJ_B,
                                                                  INCR, MAX_DATASET, N_RANGE_COEFF, PREFIX_B, i, D2_MIN_BOUND, D2_MIN_RANGE, FORECAST_SPAN)
        with open(path + "/input.txt", "w") as f_in1:
            f_in1.write(input1)
        with open(path + "/input2.txt", "w") as f_in2:
            f_in2.write(input2)
        f_out1 = open(path + "/output.txt", "w")
        f_err1 = open(path + "/error.txt", "w")
        f_out2 = open(path + "/output2.txt", "w")
        f_err2 = open(path + "/error2.txt", "w")
        f_in1 = open(path + "/input.txt", "r")
        f_in2 = open(path + "/input2.txt", "r")
        print(input1)
        print(input2)
        p1 = subprocess.Popen(prog, stderr=f_err1, stdout=f_out1, stdin=f_in1)
        p2 = subprocess.Popen(prog, stderr=f_err2, stdout=f_out2, stdin=f_in2)
        p1.wait()
        p2.wait()
        processes = []
        for j in range(0, i):
            subpath = path + "/" + str(j)
            try:
                os.mkdir(subpath)
            except:
                pass
            os.chdir(subpath)
            processes.append(Process(target=mytrainer, args=(j,i)))
            processes[j].start()
        for process in processes:
            process.join()


           
