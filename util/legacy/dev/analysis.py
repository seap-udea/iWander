from iwander import *

############################################################
#CONSTANTS AND MACROS
############################################################
#BASE DIR
BD="../"

def surrogateObjectsElements():
    wobjs=pd.read_csv(BD+"wanderer.csv")
    nominal=wobjs.iloc[0]

surrogateObjectsElements()
