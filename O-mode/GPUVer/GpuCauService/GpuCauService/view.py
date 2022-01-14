from django.http import HttpResponseRedirect
from django.http import StreamingHttpResponse,FileResponse
from django.http import HttpResponse,JsonResponse
import json
import numpy as np
import threading

import GpuCau

GpC = GpuCau.CauFullwave()


def RunGpuCau():
    GpC.Gpurun()


def GetRun(request):

    res = request.body.decode('utf-8')
    print(res)
    DataRes = json.loads(res)
    f = DataRes['f']
    nsteps = DataRes['nsteps']
    Vt = DataRes['Vt']
    Turs = np.array(DataRes['Turs'],dtype=np.float64)
    GpC._initparameters(f,nsteps,Turs,Vt)
    Thed = threading.Thread(target=RunGpuCau)
    Thed.setDaemon(True)
    Thed.start()
    return HttpResponse('ok')


def download_data(request):
    FileName = 'FullData.mat'
    FilePath = 'result/'
    file = open(FilePath+FileName, 'rb')
    response = FileResponse(file)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="{}"'.format(FileName)
    return response
def GetRunSata(request):
    res = request.GET['com']
    if res =='GetCurIterNums':
        CurNums = str(GpC.CurIterNums)
        return HttpResponse(CurNums)


