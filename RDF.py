import tools
import sys

path=sys.argv[1]
data=tools.readFile(path)
tools.RDF("./",data,50,50,15)