import os
import h5py
from pathlib import Path
import scipy.io as scio
import scipy.stats as stats
import pickle
import matplotlib.pyplot as plt
import numpy as np

class format_Stimulus_Presentation_Session_BLAES():
    def __init__(self,loc:Path,subject,plot_stimuli=False,HDF:bool=False):
        """
        Formats the 4 preprocessed filed from MATLAB into a data object
        
        Parameters
        ----------
        loc: Path
            path to preprocessed directory
        subject: str
            individual subject ID
        """
        self.root = loc.parent
        saveDir = loc/'processed'
        if not os.path.exists(saveDir):
            os.mkdir(saveDir)
        self.subject = subject
        self.name = 'stim_presentation'
        self.epoch_info = {} # intialize empty as different experiments have different requirements for this variable.

        if HDF:
            # need to write HDF5 parser. 
            # with h5py.File(loc/file, 'r') as f:
            #     # data = {key: f[key][()] for key in f.keys()}
            #     data = f['agg_signals'][()]
            data = {}
            hf = h5py.File(loc/'stimsweep.mat','r')
            keys = list(hf.keys())
            try:
                for i in keys[1:]:
                    temp = hf[i]
                    if isinstance(temp,h5py.Dataset):
                        self.seegdat = temp[0:]
                    else:
                        self.__setattr__(i,{})    
                        for k in temp.keys():
                            self.__getattribute__(i)[k] = temp[k][0:]
                temp = {}            
            # for item in self.seeg['data']:
                
            except KeyError:
                print('Mismatch in Keys')
            # for k in keys:
            #     data[k] = temp[k][0]
            # self.data = data

        else:    
            data = scio.loadmat(loc,mat_dtype=True,simplify_cells=True)
            self.data = data['signals']
        self.channels = scio.loadmat(loc/'channels.mat',mat_dtype=True,simplify_cells=True)['channels']
        self.data = {'sEEG':{k:v for k,v in zip(self.channels,self.seegdat)}}
        print('timeseries loaded')
        print(0)
        self.params = scio.loadmat(loc/'params.mat',mat_dtype=True,simplify_cells=True)['params']
        print(0)
        self.intervals = scio.loadmat(loc/'intervals.mat',mat_dtype=True,simplify_cells=True)['intervals']
        print(0)
        print(0)
        self.signalTypes = set(self.data.keys())
        print(0)
        self.fs = self.params['SamplingRate']['NumericValue']
        del self.seegdat
        # self.states = self._reshapeStates()
        self.saveSelf(fp=saveDir/'session.pkl')
        print(0)

    def saveSelf(self,fp):
        with open(fp, 'wb') as f:
            pickle.dump(self,f)
    
    @classmethod
    def class_loader(cls,filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)


    def getCommonAverages(self):
        common_avg = {}
        for sigtype in self.signalTypes:
            temp = []
            for channel,t in self.data[sigtype].items():
                temp.append(self.data[channel])
            common_avg[sigtype] = np.mean(temp,axis=0)
        return common_avg
    def _reshapeStimuliMatrix(self,stimuli):
        keys = stimuli[0].keys()
        output = {}
        for value,key in enumerate(keys):
            temp = []
            for entry in stimuli:
                temp.append(entry[key])
            temp.insert(0,{'code':value+1})
            output[key] = temp
        return output
    def _saveSelf(self,loc: Path):
        # TODO: make these modules save and load themselves from files to speed up operations, especially since this operation does not alter the data in any way.
        # https://stackoverflow.com/questions/2709800/how-to-pickle-yourself
        fname = loc/'session_aggregate.pkl'
        with open(fname, 'wb') as fp:
            pickle.dump(self)

    def epochStimulusCode_SCANtask(self,plot_states:False):
        data = self.states['StimulusCode']
        moveEpochs = {}
        onset_shift = 1000
        offset_shift = 3000
        for stim in self.stimuli.values(): # get intervals for each of the stimulus codes
            code = stim[0]['code']
            stim_type = stim[6]
            loc = np.where(data==code)
            intervals = find_intervals(loc[0])
            for i,v in enumerate(intervals):
                intervals[i] = [v[0]-onset_shift,v[1]+offset_shift]

            moveEpochs[stim_type] = intervals
        loc = np.where(data==0) # get intervals for stim code of zero (at rest)
        intervals = find_intervals(loc[0])
        for i,v in enumerate(intervals):
                intervals[i] = [offset_shift+v[0],v[1]-onset_shift]
        restEpochs = {}
        for k, int_set in moveEpochs.items():
            temp = []
            for i in int_set:
                onset = i[1]+1
                temp.append(extractInterval(intervals,onset))
            mode_int = stats.mode([i[1]-i[0] for i in temp]).mode
            for i in temp:
                if (i[1]-i[0]) != mode_int:
                    i[1] = i[0] + mode_int

            restEpochs[k] = temp                

        if plot_states:
            self.plotStimuli(moveEpochs)
        return moveEpochs, restEpochs
    
    def epochStimulusCode_screening(self,plot_states:False):
        data = self.states['StimulusCode']
        epochs = {}
        onset_shift = 0
        offset_shift = 0
        for stim in self.stimuli.values(): # get intervals for each of the stimulus codes
            code = stim[0]['code']
            stim_type = f'{stim[1]}_{code}'
            loc = np.where(data==code)
            intervals = find_intervals(loc[0])
            for i,v in enumerate(intervals):
                intervals[i] = [v[0]-onset_shift,v[1]+offset_shift]
            epochs[stim_type] = intervals
        if plot_states:
            self.plotStimuli(epochs)
        return epochs
    def _reshapeStates(self):
        states = {}
        for state, val in self.states.items():
            if np.all(val[0]==np.max(val[0])):
                pass # Exclude states that do not contain any information, ie array contains all of one value. 
            else:
                states[state] = val[0]
        return states
    def plotStimuli(self,epochs):
        import distinctipy
        data = self.states['StimulusCode']
        t = np.linspace(0,len(data)/2000,len(data))
        fig=plt.figure()
        ax = plt.subplot(1,1,1)
        ax.plot(t,data)
        cmap = distinctipy.get_colors(len(epochs))
        for i,k in enumerate(epochs.keys()):
            val = int(k.split('_')[-1])
            for q,j in enumerate(epochs[k]):
                if q == 0:
                    # ax.axvline(t[j[0]], c=cmap[i], label=k)
                    # ax.axvline(t[j[1]], c=(0,0,0),label='_')
                    ax.plot(t[j[0]:j[1]],val*np.ones(len(t[j[0]:j[1]])),c=cmap[i], label=k)
                else:
                    # ax.axvline(t[j[0]], c=cmap[i],label='_')
                    # ax.axvline(t[j[1]], c=(0,0,0),label='_')
                    ax.plot(t[j[0]:j[1]],val*np.ones(len(t[j[0]:j[1]])),c=cmap[i], label="_")
        # scio.savemat(self.root/'analyzed'/'stimcode.mat',{'stim':data})
        # scio.savemat(self.root/'analyzed'/'stimuli.mat',epochs)
        ax.set_ylim(-1,max(data)+2)
        ax.legend()
        ax.set_title(self.name)
        plt.show()
        
def stimulus_presentation_loader(fp,**kwargs)->format_Stimulus_Presentation_Session_BLAES:
    objPath = fp/'processed'/'session.pkl'
    if os.path.isfile(objPath):
        return format_Stimulus_Presentation_Session_BLAES.class_loader(objPath)
    else:
        return format_Stimulus_Presentation_Session_BLAES(**kwargs)
def find_intervals(array):
    # Initialize the list of intervals
    intervals = []
    # Start the current interval with the first element
    start = array[0]

    # Iterate over the array starting from the second element
    for i in range(1, len(array)):
        # If the difference from the previous to the current is not 1, we've found the end of an interval
        if array[i] - array[i - 1] != 1:
            # Add the (start, end) of the interval to the list
            intervals.append((start, array[i - 1]))
            # Start a new interval with the current element
            start = array[i]
    
    # Add the last interval if the array does not end on a jump based on length of previous intervals
    if len(intervals) > 0:
        int_length = intervals[-1][1] - intervals[-1][0]
        intervals.append((start, start+int_length))
    else:
        intervals = [(array[0],array[-1])]
    return intervals


class screening_session(object):
    """docstring for screening_session."""
    # TODO: make these modules save and load themselves from files to speed up operations, especially since this operation does not alter the data in any way.
    # https://stackoverflow.com/questions/2709800/how-to-pickle-yourself
        
    def __init__(self, loc:Path,subject:str,plot_stimuli: bool=False, HDF:bool=True):
        super(screening_session, self).__init__()
        self.motor =        format_Stimulus_Presentation_Session_BLAES(loc/'motor/preprocessed',subject,plot_stimuli=plot_stimuli,HDF=HDF)
        self.motor.name = 'motor'
        self.sensory =      format_Stimulus_Presentation_Session_BLAES(loc/'sensory/preprocessed',subject,plot_stimuli=plot_stimuli,HDF=HDF)
        self.sensory.name = 'sensory'
        self.sensorimotor = format_Stimulus_Presentation_Session_BLAES(loc/'sensory-motor/preprocessed',subject,plot_stimuli=plot_stimuli,HDF=HDF)
        self.sensorimotor.name = 'sensorimotor'

def extractInterval(intervals,b):
    for i in intervals:
        if i[0] ==b: 
            return i
    return None

def sliceArray(array, interval):
    return array[interval[0]:interval[1]]

if __name__ == '__main__':
    loc = Path(r'/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param/BJH050')
    subject = 'BJH050'
    kwargs = {'loc':loc,'subject':subject,'plot_stimuli':False,'HDF':False}
    a = stimulus_presentation_loader(loc,loc=loc,subject=subject,plot_stimuli=False,HDF=True)
    a.states = a._reshapeStates()
    breakpoint()
    print(0)
