import sys
import os

def generatePocosInput(filename, args):
    hexaedra, tetrahedra, nx, ny, nz = args

    with open(filename, "w") as file:
        fileString = f'''{hexaedra}
{tetrahedra}
{nx} {ny} {nz}
1 1 1:{nz}
{nx} 1 1:{nz}
1 {ny} 1:{nz}
{nx} {ny} 1:{nz}'''

        file.write(fileString)
if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print("Incorrect number of arguments")
        os.exit(1)
    
    generatePocosInput(sys.argv[1:])






