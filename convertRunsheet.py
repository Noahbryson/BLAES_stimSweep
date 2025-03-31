import pandas as pd
import openpyxl

inputPath = r'/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param/ParamSweep_ALL_PATIENTS_DescriptiveInfo.csv'
outputPath = r'/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param/Subject_Locations.xlsx'
with open(inputPath,'r') as fp:
      df = pd.read_csv(fp)
      # workbook = openpyxl.load_workbook(fp)
      # sheet=workbook.active



def sparseStimInformation(ID:str,stims:str):
      if ID.find('BJH')>=0:
            site='BJH'
      elif ID.find('UIC')>=0:
            site='UIC'
      else:
            site='NAN'
      pairs = stims.split(' and ')
      pair1 = pairs[0].replace(' - ','-')
      pair2 = pairs[1].replace(' - ','-')
      pair1 = pair1.replace('-','_')
      pair2 = pair2.replace('-','_')
      pair1 = pair1.replace('(_)','(-)')
      pair2 = pair2.replace('(_)','(-)')
      
      triggers = '_'.join([pair1,pair2])
      triggers = ','.join(triggers.split('_'))
      triggers = triggers.replace('(-)','')
      triggers = triggers.replace('(+)','')   
      return pair1,pair2,triggers,site
def addFillerCols(df:pd.DataFrame):
      fillers = ['' for i in range(len(df))]
      cols=['loc1','side1','loc2','side2']
      for i in cols:
            df[i] = fillers
      return df
def segmentDf(df,targets):
      output=pd.DataFrame()
      for i in targets:
            output[i]=df[i]
      return output


df_col_order = ['Subject','Pair1','Pair2','Triggers','site']
df = df.dropna(subset='PATIENT ID')
df[['Pair1','Pair2','Triggers','site']] = df.apply(lambda x: sparseStimInformation(x['PATIENT ID'],x['AMYGDALA STIM CONTACTS']),axis=1,result_type='expand')
# df = addFillerCols(df)
df.rename({'PATIENT ID':'Subject'},axis=1,inplace=True)
output = segmentDf(df,df_col_order)
output.set_index('Subject',inplace=True)
output.to_excel(outputPath)

print(0)