import numpy as np
from scipy.spatial import distance

'''
Requires the binarized signal as an array of shape (number of frames, number of cells)
and cellular positions (x,y) as an array of shape (number of cells, 2)

A folder with the name 'waves' must be present in the folder where this script is located.
'''

#Loads required data
binsig=np.loadtxt("Fbinsig.txt") #binarized cell signals. 1: during Ca2+ oscillations, otherwise 0
                                #signal is of shape: (frame index, cell index)
pos=np.loadtxt("Fkoordinate.txt") #cell coordinates. Array of shape: (number of cells, 2) (x,y) columns

##Set required parameters
sampling=10.0 #sampling rate during Ca2+ signal acquisition in Hz
Tth=0.7 
'''Tth -> thershold time (seconds): if the activation times of two cells is lower than Tth
    the signal is considered to have jumped from one cell to the other 
    (if the cells are within a preset threshold distance). The threshold distance (Rth) is calculated
    based on the average intercellular distance and a fraction of the std: line 34
'''

##############Dont change things from here on########################

act_sig=np.zeros(binsig.shape, int)
event_num=[]
NC=len(pos) #stevilo celic
frameTh=int(round(Tth*sampling))

R=distance.cdist(pos, pos, 'euclidean')
Ravg = np.average(R)
Rstd = np.std(R)
Rth=Ravg - 0.15*Rstd
##########################################################
##################WAVE DETECTION###########################


        

print(f'Average intercellular distance is {Ravg} with a standard deviation of {Rstd}')
print(f'Set threshold distance for signal propagation is {Rth}')

###Finds all neighbours of cells which are within the threshold distance
neighbours=[]
for j in range(len(pos)):
    neighbours.append(np.where((R[:,j]<=Rth) & (R[:,j]!=0))[0])

##Finds all frames which contain non-zero values (ones) <--- binarized activity  
nonzero={}
frame=np.where(binsig==1)[0] #poisce vse frejme kjer je kakšna celica aktivna
for v in frame: #zanka po frejmih z aktivnostjo
    nonzero[v]=list(np.where(binsig[v,:]==1)[0]) #v frejmu z aktivnostjo poisce vse celice, ki so dejansko aktivnee

counter=0
for i in nonzero.keys(): #loop over all frames with cellular activity
    if i%50==0: print(i)
    if(counter==0):  #analyses first frame
        k=1
        for j in nonzero[i]: #loop over active cells in this frame
            act_sig[i,j]=k
            k+=1
        
        for nn in nonzero[i]:
            for nnn in list(set(neighbours[nn]).intersection(set(nonzero[i]))):
                act_sig[i,nn]=min(act_sig[i,nn],act_sig[i,nnn])
                act_sig[i,nnn]=act_sig[i,nn]
        
        un_num=np.unique(act_sig[i,:])      
        event_num=list(set(event_num).union(set(un_num)))
                
        max_event_num=max(event_num)
        counter+=1
        
    else: ##analyses all the other frames with activity
        k=max_event_num+1
        
        for j in nonzero[i]: #loop over active cells in frame
            if(binsig[i-1,j]==0): #if the cell is active in this frame and not in the previous it assignes a new activity index
                act_sig[i,j]=k
                k+=1
            else: #else it assignes previous activity index
                act_sig[i,j]=act_sig[i-1,j]
        
        #loops over all cells and neighbours and clusters cells which are close in space (Rth) and time (Tth)
        for nn in nonzero[i]:
            for nnn in list(set(neighbours[nn]).intersection(set(nonzero[i]))):
                if(act_sig[i,nn]!=0 and act_sig[i,nnn]!=0 and act_sig[i-1,nn]!=0 and np.sum(binsig[i-frameTh:i+1,nn])<=frameTh and act_sig[i-1,nnn]==0 and nn!=nnn):
                    act_sig[i,nnn]=act_sig[i,nn]
                elif(act_sig[i,nn]!=0 and act_sig[i,nnn]!=0 and act_sig[i-1,nn]==0 and act_sig[i-1,nnn]!=0 and np.sum(binsig[i-frameTh:i+1,nnn])<=frameTh and nn!=nnn):
                    act_sig[i,nn]=act_sig[i,nnn]
                elif(act_sig[i,nn]!=0 and act_sig[i,nnn]!=0 and act_sig[i-1,nn]==0 and act_sig[i-1,nnn]==0 and nn!=nnn):
                    act_sig[i,nn]=min(act_sig[i,nn],act_sig[i,nnn])
                    act_sig[i,nnn]=act_sig[i,nn]
        
                    
        un_num=np.unique(act_sig[i,:])
        event_num=list(set(event_num).union(set(un_num)))
                
        max_event_num=max(event_num)
        counter+=1

##saves the data as 'act_sig.txt' in the folder 'waves'  
'''
The data in act_sig.txt is of the same shape as binsig (binarized cellular activity)
shape: (number of frames, number of cells)

The values in the array are integers (>0) and correspond to the event (event number) 
the cell is considered part of. 
Example: all cells which have the number 54 belong to the event/wave with index
54. This means that you can extract all cells and their corresponding frames with
this index and reconstruct the wave propagation.
'''
np.savetxt("waves//act_sig.txt",act_sig,fmt='%d')



    
    











