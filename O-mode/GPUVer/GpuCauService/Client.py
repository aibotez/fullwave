import requests,time,os,json

url = 'http://222.204.54.82:6000/'
def initparameters(f,nsteps,Vt,Turs,url0):
    data = {
        'f':f,
        'nsteps':nsteps,
        'Vt':Vt,
        'Turs':Turs
    }
    parameters = json.dumps(data)
    url = url0+'run/'
    res = requests.post(url,data=parameters)
    if res.text == 'ok':
        print('命令发起成功')
def GetRunSata(url0):
    url = url0+'GetRunSata/?com=GetCurIterNums'
    sata = ''
    while 'finsh' not in sata:
        sata1 = requests.get(url).text
        if sata1 != sata:
            sata = sata1
            print(sata)
        time.sleep(1)
def GetResultData(url0):
    curtrace = os.getcwd()
    FilePath = 'GetResult/'
    if os.path.isdir(curtrace + '/{}'.format(FilePath)):
        pass
    else:
        os.makedirs(curtrace + '/{}'.format(FilePath))

    url = url0+'download_data/'
    r = requests.get(url)
    FileName = 'FullData.mat'
    with open(FilePath+FileName, "wb") as code:
        code.write(r.content)






f = 30*10**9
c = 3*10**8

nsteps=30000
Vt = 1 * c / 20
Turs = []


initparameters(f,nsteps,Vt,Turs,url)
GetRunSata(url)
GetResultData(url)