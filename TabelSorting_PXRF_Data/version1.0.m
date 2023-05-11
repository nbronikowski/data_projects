% # First in bash:
% # grep -n File org.csv | awk -F: '{print $1}' > lines
% #
% CREATES A FILE WITH FILE NUMBERS FOR INDEX 

% # Then in python:
% import pandas as pd
% lines = pd.read_csv('lines.csv')
% dataList = []
% for i in range(len(lines)-1):
%     iHeaderLine = lines.iloc[i][0]
%     iNextHeader = lines.iloc[i+1][0]
%     print(i, iHeaderLine,iNextHeader-iHeaderLine-1)
%     dataList.append(pd.read_csv('org.csv', dtype=str, skiprows=iHeaderLine-1,nrows=iNextHeader-iHeaderLine-1))
% data = pd.concat(dataList)
% data.to_csv('final.csv') (edited) 



