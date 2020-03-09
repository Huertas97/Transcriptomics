import os

def StatesFileCreator(Directory):
    os.chdir(Directory)
    Files = os.listdir()
    for File in Files:
        with open(File, "r") as file:
            Output1 = open("States_E1.bed", "a")
            Output2 = open("States_E2.bed", "a")
            istart = 0
            iend = 200
            for line in file:
                probs = []
                if "chr" in line:
                    chromosome = line.strip().split("\t")[1]
                elif "E" not in line:
                    probs = line.strip().split("\t")
                    maximum = 0.0
                    for pos, element in enumerate(probs):
                        prob = float(element)
                        if prob > maximum:
                            maximum = prob
                            state = "E" + str(pos+1)
                    if maximum >= 0.9 and state == "E1":
                        Output1.write(chromosome + "\t" + str(istart) + "\t" + str(iend) + "\t" + state + "\n")
                    if maximum >= 0.9 and state == "E2":
                        Output2.write(chromosome + "\t" + str(istart) + "\t" + str(iend) + "\t" + state + "\n")
                    istart += 200
                    iend += 200
    Output1.close()
    Output2.close()
    os.chdir("../")

StatesFileCreator("./Monocyte1")
StatesFileCreator("./Monocyte2")
