import requests
import re
import numpy as np
import pandas as pd
import json


def getnist(URL_GET):
    # URL_GET = "https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=50&Formula=&gtype=4&range=U&lower=&upper=&density=&frames=no"
    #构建请求参数
    params = {}
    #发送请求
    response = requests.get(URL_GET,params=params)
    # copy?
    t = response.text

    # inihead
    head = {}

    # element name
    b1 = t.find('<b>')
    b_name = t[b1:].find('&')
    head['name'] = t[b1+8:b1+b_name]

    # Z
    b_Z = t[b1:].find('Z =') + b1
    b_Zb = t[b_Z:].find(')')
    head['Z'] = int(t[b_Z+3:b_Z+b_Zb])

    # weight
    b1 = t[b1:].find('Atomic weight:') + b1
    b_w = t[b1:].find('=') + b1
    b_wb = t[b1:].find('g mol') + b1
    head['weight'] = float(t[b_w+4:b_wb-1])

    # density
    b1 = t[b1:].find('Nominal density:') + b1
    b_d = t[b1:].find('=') + b1
    b_db = t[b1:].find('\n') + b1
    head['density'] = float(t[b_d+2:b_db])

    # edges
    b_e = t[b1:].find('edge')
    b1 = b_e + b1
    head['edges'] = int(t[b1 - 4:b1 - 1])

    # Edge energies
    b1 = t[b1:].find('Edge energies:') + b1
    b_e = t[b1:].find('<pre>') + b1
    b_ed = t[b1:].find('</pre>') + b1
    head['E_edges'] = [float(x) for x in re.findall('\d+.\d+E[+-]\d+', t[b_e:b_ed])]

    # relativistic correction
    b1 = t[b1:].find('Relativistic correction estimate') + b1
    b_r1 = t[b1:].find('=') + b1
    b_r1d = t[b_r1:].find(',') + b_r1
    head['f_relH82'] = float(t[b_r1+1:b_r1d])
    b_r2 = b_r1d
    b_r2d = t[b_r2:].find('<i>') + b_r2
    head['f_rel35CL'] = float(t[b_r2+1:b_r2d])

    # Nuclear Thomson correction
    b1 = t[b1:].find('Nuclear Thomson correction') + b1
    b_nt = t[b1:].find('=') + b1
    b_ntd = t[b_nt:].find('<i>') + b_nt
    head['f_NT'] = float(t[b_nt+1:b_ntd])

    # data table
    b1 = t[b1:].find('Form Factors') + b1
    data = [float(x) for x in re.findall('\d+.\d+E[+-]\d+', t[b1:])]
    data = np.array(data).reshape(-1, 8)
    data_df = pd.DataFrame(data, columns=['E', 'f1', 'f2', 'PE', 'C+i', 'Total', 'K', 'lambda'])

    # return
    return head, data_df

path = './elements/'
for Z in range(1,92):
    URL_GET = "https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=" + str(Z) + \
              "&Formula=&gtype=4&range=U&lower=&upper=&density=&frames=no"
    [head, data] = getnist(URL_GET)
    name = head['name']
    with open(path+name+'.txt', 'w') as f:
        json.dump(head, f, indent=3)
    data.to_csv(path+name+'.csv', index=None)