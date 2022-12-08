import re

def extractData2(file):
	numberRegex = re.compile(r'-?\d+.\d+')
	coordStringValues = []
	zcornStringValues = []
	numberOfCellsPerAxis = []
	for linha in file:
		if linha.lstrip().startswith("NX"):
			linha = next(file)
			numberOfCellsPerAxis.append(linha)
		if linha.startswith("COORD"):
			while not linha.startswith("ZCORN"):
				coordStringValues.append(linha)
				linha = next(file)
			while not linha.startswith("DTOP"):
				zcornStringValues.append(linha)
				#print("read: " + linha)
				linha = next(file)
	coordValues = [float(values) for line in coordStringValues for values in numberRegex.findall(line)]
	zcornValues = [float(values) for line in zcornStringValues for values in numberRegex.findall(line)]
	numberOfCellsPerAxisValues = [float(values) for line in numberOfCellsPerAxis for values in re.compile(r'\d+').findall(line)]

	return coordValues, zcornValues, numberOfCellsPerAxisValues

def extractData3(file):
	numberRegex = re.compile(r'-?\d+.\d+')
	end_of_file = False
	coordStringValues = []
	zcornStringValues = []
	numberOfCellsPerAxis = []
	for linha in file:
		if linha.lstrip().startswith("NX"):
			linha = next(file)
			numberOfCellsPerAxis.append(linha)
		if linha.lstrip().upper().startswith("COORD"):
			while not linha.lstrip().upper().startswith("ZCORN"):
				coordStringValues.append(linha)
				linha = next(file)
			while not linha.lstrip().upper().startswith("PORO"):
				zcornStringValues.append(linha)
				try:
					linha = next(file)
				except StopIteration:
					break
			break
	coordValues = extractCoord(coordStringValues)
	zcornValues = extractCoord(zcornStringValues)
	numberOfCellsPerAxisValues = [int(values) for line in numberOfCellsPerAxis for values in re.compile(r'\d+').findall(line)]

	return coordValues, zcornValues, numberOfCellsPerAxisValues

def getInactiveCells(file):
	mult = re.compile(r'(\d+)\*(\d+)')

	binary_list = []
	with open(file) as f:
		for linha in f:
			linha = linha.strip()
			if linha.upper().startswith("NULL"):
				continue
			split_linha = linha
			list_tokens = split_linha.split(" ")
			tokenize = []
			for token in list_tokens:
				mo = mult.search(token)
				if mo != None:
					for i in int(mo.group(1)) * mo.group(2):
						tokenize.append(int(i))
				else:
					tokenize.append(int(token))
					
			binary_list.extend(tokenize)
	
	return binary_list

def getPermeability(file):
	mult = re.compile(r'(\d+)\*(\d+.\d+)')

	permeability_list = []
	with open(file) as f:
		for linha in f:
			linha = linha.strip()
			if linha.upper().startswith("NULL"):
				continue
			split_linha = linha
			list_tokens = split_linha.split(" ")
			tokenize = []
			for token in list_tokens:
				mo = mult.search(token)
				if mo != None:
					for i in range(int(mo.group(1))):
						tokenize.append(float(mo.group(2)))
				else:
					try:
						tokenize.append(float(token))
					except:
						pass
						
					
			permeability_list.extend(tokenize)
	
	return permeability_list

def extractCoord(file):
	mult = re.compile(r'(\d+)\*(-?\d+(.\d+)?)')
	coordValues = []
	for linha in file:
		if linha.lstrip().upper().startswith("COORD") or linha.lstrip().upper().startswith("ZCORN"):
			continue
		
		list_tokens = linha.lstrip().split("\t")
		
		if(len(list_tokens) == 1): # If tabs don't work split with spaces 
			list_tokens = list_tokens[0].split(" ")
		
		for token in list_tokens:
			mo = mult.search(token)
			if mo != None:
				for i in range(int(mo.group(1))):
					coordValues.append(float(mo.group(2)))
			else:
				if token != "\n" and token != "" and token != "/\n":
					coordValues.append(float(token))
	return coordValues

def extractHeader(file):
	#grid_200x200x49.asc
	numberCells = []
	fileNames = []
	name = re.compile(r'\'(.*)\'')
	for linha in file:
		if linha.lstrip().startswith("SPECGRID"):
			linha = next(file)
			numberCells = linha.lstrip().split(" ")
		elif linha.lstrip().startswith("INCLUDE"):
			linha = next(file)
			fileNames.append(name.search(linha).group(1))

	return [int(numberCells[0]), int(numberCells[1]), int(numberCells[2])], fileNames

def extractDataDTOP(file):
	#numberRegex = re.compile(r'-?\d+\.\d+')
	numberRegex2 = re.compile(r'(\d+(\.\d+)?)')
	numberOfCellsPerAxis = []
	dtopValues = []
	width = []
	for linha in file:
		if linha.lstrip().startswith("NX"):
			linha = next(file)
			numberOfCellsPerAxis.append(linha)
		elif linha.lstrip().startswith("DX"):
			linha = next(file)
			width.append(linha)
		elif linha.lstrip().startswith("DTOP"):
			while True:
				try:
					linha = next(file)
					dtopValues.append(linha)
				except Exception as e:
					 break

	# print(numberRegex2.findall(width[0]))
	numberOfCellsPerAxisValues = [int(values) for line in numberOfCellsPerAxis for values in re.compile(r'\d+').findall(line)]
	#dtopValues = [float(values[0]) for line in dtopValues for values in numberRegex2.findall(line)]
	dtopValues = extractCoord(dtopValues)
	width = [float(values[0]) for line in width for values in numberRegex2.findall(line)]
	return numberOfCellsPerAxisValues, dtopValues, width


