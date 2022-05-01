A simple algorithm for Ca2+ wave detection in acute pancreatic tissue slices
based on the binarized Ca2+ activity of beta cells in Islets of Langerhans
recorded with high spatio-temporal resolution confocal laser scaning microscopy.

Example binarized beta cell activity data is in the 'Fbinsig.txt' file with corresponding
cellular coordinates in the 'Fkoordinate.txt' file. 

Fbinsig.txt --> array like structure with shape (number of recorded frames, number of cells)
Fkoordinate.txt --> array like structure with shape (number of cells, 2): x, y coordinates 

Algorithm is in the file 'wave_detection.py' and takes the beforementioned .txt files
+ a few parameters which are set in the file (explanations are in the file).

Result is of the same shape as the binarized activity array with values corresponding to the
waves each cell belongs to. Example: all cells with value 52 in the resulting array belong to
the same event with number 52. Each wave can therefore be reconstructed by extracting cells 
that belong to the same event at a given time.