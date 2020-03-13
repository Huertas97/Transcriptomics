import os
import sys

def StatesFileCreator(Directory, Output1Name, threshold):
    if len(sys.argv) < 4:
        print("Es obligatorio pasar los siguientes argumentos en el orden especificado: nombre del directorio de inicio, nombre del archivo de salida",
              "correspondiente al estado E1 y E2 (sin el _E1.bed o _E2.bed), threshold de probabilidad.")
        exit()
    threshold = float(threshold)
    os.chdir(Directory)
    Files = os.listdir()
    for File in Files:
        with open(File, "r") as file:
            Output1 = open(str(Output1Name)+"_E1.bed", "a")
            Output2 = open(str(Output1Name)+"_E2.bed", "a")
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
                    if maximum >= threshold and state == "E1":
                        Output1.write(chromosome + "\t" + str(istart) + "\t" + str(iend) + "\t" + state + "\t" + str(maximum) + "\n")
                    if maximum >= threshold and state == "E2":
                        Output2.write(chromosome + "\t" + str(istart) + "\t" + str(iend) + "\t" + state + "\t" + str(maximum) + "\n")
                    istart += 200
                    iend += 200
    print("Creando el fichero", str(Output1Name)+"_E1.bed", "usando como threshold la probabilidad de estado =", threshold)
    print("Creando el fichero", str(Output1Name)+"_E2.bed", "usando como threshold la probabilidad de estado =", threshold)
    Output1.close()
    Output2.close()
    os.chdir("../")
    print()

StatesFileCreator(Directory = sys.argv[1], Output1Name = sys.argv[2], threshold = sys.argv[3])
