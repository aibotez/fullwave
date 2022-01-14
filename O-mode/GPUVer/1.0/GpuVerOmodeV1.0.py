from numba import cuda # 从numba调用cuda
cuda.select_device(0)
import numpy as np
import math,os
import matplotlib.pyplot as plt
import time,random
import scipy.io as sio
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32

@cuda.jit
def gpu_print(N):
    idx = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x
    if (idx < N):
        print(idx)
def OneDtoTwoD(n,ny):
    j = int(n / ny)
    i = n % ny
    return [i, j]

@cuda.jit
def Dne(omega_pe,ne0,vt,ntp,dt,ddx,ddy,nx,ny,Turs):
    i = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    j = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if i < nx and j < ny:
        t = (ntp - 1) * dt
        yc = 0.32
        yw = 0.1#10*lambda
        Atur = 0
        e = 1.602 * 10 ** (-19)
        me = 9.1096 * 10 ** (-31)
        epsilon0 = 8.85 * 10 ** (-12)
        nesm = 0
        xi = i * ddx
        at = ne0[i, j] * 1 / 100
        dyw = math.pow((j*ddy-yc)/yw,2)
        for kxi in Turs:
            Atur = at * math.exp(-dyw)*math.sin(kxi * (xi - vt * t))
            nesm = nesm + Atur
        ls = (ne0[i, j] + nesm) * e ** 2 / epsilon0 / me
        omega_pe[i, j] = ls
@cuda.jit
def Dne0(omega_pe,ne0,kx,vt,ntp,dt,ddx,ddy,nx,ny):
    i = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    j = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if i < nx:
        if j < ny:
            t = (ntp-1) * dt
            Atur = 0
            e = 1.602 * 10**(-19)
            me = 9.1096 * 10**(-31)
            epsilon0 = 8.85 * 10**(-12)
            nesm = 0
            xi = i * ddx
            at = ne0[i,j] * 1 / 100
            kxi = 1 * kx / 3
            Atur = at * math.sin(kxi * (xi - vt * t))
            nesm = nesm + Atur
            kxi = 2 * kx / 3
            Atur = at * math.sin(kxi * (xi - vt * t))
            nesm = nesm + Atur
            kxi = 3 * kx / 3
            Atur = at * math.sin(kxi * (xi - vt * t))
            nesm = nesm + Atur
            ls = (ne0[i,j] + nesm) * e**2 / epsilon0 / me
            omega_pe[i,j] = ls

@cuda.jit
def setAntennaBoundary(AntennaPoint,rowAn,dzG):
    i = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    j = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if i <rowAn:
        dzG[AntennaPoint[i,0],AntennaPoint[i,1]] = 0
@cuda.jit
def AddSource(dzG,lenaddS,idx0,idx3,ddx,xatp,f,dt,t0,T,w,Launchdata):
    i = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    j = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if (i<lenaddS):
        UA=5
        Uofft=100
        n = T-1
        Ui = UA * (math.exp(n / Uofft) - math.exp(-n / Uofft)) / (math.exp(n / Uofft) + math.exp(-n / Uofft))
        idxi = idx0 + i
        idxj = idx3
        pl = idxi * ddx - xatp
        signal = Ui * math.exp(-pl ** 2 / (w ** 2)) * math.cos(-2 * np.pi * f * dt * (t0 - T))
        dzG[idxi-1,idxj-1] = dzG[idxi-1,idxj-1]+signal
        if i==0 and j==0:
            Launchdata[T-1] = Ui * math.cos(-2 * np.pi * f * dt * (t0 - T))

@cuda.jit
def GetReceiveData(ezG,idx1,idx2,idx3,Receives,n):
    i = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    j = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if i < idx2-idx1+1:
        idxi = idx1 + i
        idxj = idx3
        Receives[i,n-1] = ezG[idxi-1,idxj-1]


@cuda.jit
def Caudz(dzG, hyG, hxG,nx,ny,gi3,gj3,gi2,gj2,ddx,ddy):
    ii = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    jj = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if ii>=1 and jj>=1:
        if ii < nx and jj < ny:
            dzG[ii,jj] = gi3[ii] * gj3[jj] * dzG[ii,jj] + 0.5 * gi2[ii] * gj2[jj] * (
                    hyG[ii,jj] - hyG[ii-1,jj] - (ddx / ddy) * ((hxG[ii,jj] - hxG[ii,jj-1])))
@cuda.jit
def Cauez(ezG,dzG,sxG,sx1G,sx2G,omega_peG,dt,vc,nx,ny):
    epd = math.exp(-vc * dt)
    ii = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    jj = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if ii >= 1 and jj >= 1:
        if ii < nx and jj < ny:
            ezG[ii,jj] = dzG[ii,jj] - sxG[ii,jj]
            sxG[ii,jj] = (1 + epd) * sx1G[ii,jj] - epd * sx2G[ii,jj] + (omega_peG[ii,jj] * dt / vc) * (1 - epd) * ezG[ii,jj]
            sx2G[ii,jj] = sx1G[ii,jj]
            sx1G[ii,jj] = sxG[ii,jj]
@cuda.jit
def CauHxy(ez,ihx,ihy,hx,hy,fj3,fj2,fi1,fi3,fi2,fj1,ddx,ddy,nx,ny):
    ii = cuda.blockDim.x * cuda.blockIdx.x + cuda.threadIdx.x
    jj = cuda.blockDim.y * cuda.blockIdx.y + cuda.threadIdx.y
    if ii < nx-1 and jj < ny-1:
        curl_e = ez[ii,jj] - ez[ii,jj+1]
        ihx[ii,jj] = ihx[ii,jj] + curl_e
        hx[ii,jj] = fj3[jj] * (ddx / ddy) * hx[ii,jj] + fj2[jj] * 0.5 * (ddx / ddy) * (curl_e + fi1[ii] * ihx[ii,jj])
        curl_e = ez[ii,jj] - ez[ii+1,jj]
        ihy[ii,jj] = ihy[ii,jj] + curl_e
        hy[ii,jj] = fi3[ii] * hy[ii,jj] - fi2[ii] * (0.5 * curl_e + fj1[jj] * ihy[ii,jj])

class CauFullwave():
    def __init__(self):
        self.c = 3*10^8
        self.ddx=0
        self.ddy=0
        self.vc =0
        self.dt =0
        self.Vt = None
        self.nx=0
        self.ny=0
        self.dzG=None
        self.ezG = None
        self.hxG = None
        self.hyG = None
        self.sxG = None
        self.sx1G = None
        self.sx2G = None
        self.fj3G=None
        self.fi3G=None
        self.fi2G = None
        self.fi1G = None
        self.fj2G = None
        self.fj1G = None
        self.gj2G = None
        self.gj3G = None
        self.gi2G = None
        self.gi3G = None
        self.omega_peG = None
        self.threads_per_block=None
        self.blocks_per_grid=None
        self.c = 3*10**8
        self.f = None
        self.R = None
        self.Z = None
        self.lambda0 = None
        self.k0 = None
        self.nsteps = None
        self.ne = None
        self.neG=None
        self.xatp = 0.2
        self.Ax = 0.245
        self.AtL = 0.05
        self.WGw = 0.01 / 2
        self.AntennaPoints = None
        self.AntennaPointsG=None
        self.rowAn=None
        self.LaunchData = None
        self.ReceiveDatas = None
        self.ReceiveData = None
        self.LaunchDataG = None
        self.ReceiveDatasG = None
        self.t = None

        self.phase = None
        self.tphase = None
        self.ReflI = None
        self.ReflQ = None
        self.TurGs = None


    def _initparameters(self,f,nsteps,Turs,Vp):
        print('....../初始化参数/.......')
        self.lambda0 = self.c/f
        self.f = f
        self.Vt = Vp
        self.ddx = self.lambda0/20
        self.ddy = self.lambda0/20
        self.R = np.linspace(0,0.4,int(0.4/self.ddx)+1).astype(np.float64)
        self.Z = np.linspace(0, 0.5, int(0.5/self.ddy)+1).astype(np.float64)
        self.vc = 5*10**(-3)
        self.nx=np.size(self.R)
        self.ny=np.size(self.Z)
        nx = self.nx
        ny = self.ny
        self.dt = self.ddx / self.c / 2
        self.LaunchData = np.zeros(nsteps).astype(np.float64)

        self.nsteps = nsteps
        self.ne = self.Caune(self.R,self.Z)
        LaunchAntennaPoint=self.setAntenna(self.xatp, 1.5 * self.AtL,self.WGw, 1.3 * self.AtL, self.ddx, self.ddy)
        ReceiveAntennaPoint = self.setAntenna(self.Ax, 1.5 * self.AtL,self.WGw, 1.3 * self.AtL, self.ddx, self.ddy)
        [nxLaunch,ignor] = LaunchAntennaPoint.shape
        [nxReceive, ignor] = ReceiveAntennaPoint.shape
        self.AntennaPoints = np.zeros((nxLaunch+nxReceive,2)).astype(np.int32)
        self.AntennaPoints[0:nxLaunch,:] = LaunchAntennaPoint
        self.AntennaPoints[nxLaunch:nxLaunch+nxReceive, :] = ReceiveAntennaPoint
        self.rowAn = int(nxLaunch+nxReceive)
        self.omega_pe = self.omega_pecal(self.ne)

        dz = np.zeros((nx,ny)).astype(np.float64)
        ez = np.zeros((nx,ny)).astype(np.float64)
        ihx = np.zeros((nx,ny)).astype(np.float64)
        ihy = np.zeros((nx,ny)).astype(np.float64)
        hx = np.zeros((nx,ny)).astype(np.float64)
        hy = np.zeros((nx,ny)).astype(np.float64)
        sx = np.zeros((nx,ny)).astype(np.float64)
        sx1 = np.zeros((nx,ny)).astype(np.float64)
        sx2 = np.zeros((nx,ny)).astype(np.float64)
        # omega_pe = np.zeros((nx,ny)).astype(np.float32)
        #
        gi3 = np.ones(nx).astype(np.float64)
        gj3 = np.ones(ny).astype(np.float64)
        gi2 = np.ones(nx).astype(np.float64)
        gj2 = np.ones(ny).astype(np.float64)
        fj3 = np.ones(ny).astype(np.float64)
        fj2 = np.ones(ny).astype(np.float64)
        fi1 = np.zeros(nx).astype(np.float64)
        # fi3 = np.zeros(nx)
        # fi2 = np.zeros(nx)
        # fj1 = np.zeros(nx)
        fi3 = np.ones(nx).astype(np.float64)
        fi2 = np.ones(nx).astype(np.float64)
        fj1 = np.zeros(ny).astype(np.float64)
        coeff_pml = 0.33
        npml = 8
        iternpml =np.linspace(1,8,8).astype(np.int32)
        for i in iternpml:
            xnum = npml - i
            xd = npml
            xxn = xnum / xd
            xn = coeff_pml * xxn**3
            gi2[i - 1] = 1 / (1 + xn)
            gi2[nx - 1 - i - 1] = 1 / (1 + xn)
            gi3[i - 1] = (1 - xn) / (1 + xn)
            gi3[nx - 1 - i - 1] = (1 - xn) / (1 + xn)
            gj2[i - 1] = 1 / (1 + xn)
            gj2[ny - 1 - i - 1] = 1 / (1 + xn)
            gj3[i - 1] = (1 - xn) / (1 + xn)
            gj3[ny - 1 - i - 1] = (1 - xn) / (1 + xn)
            xxn = (xnum - 0.5) / xd
            xn = coeff_pml * xxn**3
            fi1[i - 1] = xn / 2
            fi1[nx - 2 - i - 1] = xn / 2
            fi2[i - 1] = 1 / (1 + xn)
            fi2[nx - 2 - i - 1] = 1 / (1 + xn)
            fi3[i - 1] = (1 - xn) / (1 + xn)
            fi2[nx - 2 - i - 1] = (1 - xn) / (1 + xn)
            fj1[i - 1] = xn / 2
            fj1[ny - 2 - i - 1] = xn / 2
            fj2[i - 1] = 1 / (1 + xn)
            fj2[ny - 2 - i - 1] = 1 / (1 + xn)
            fj3[i - 1] = (1 - xn) / (1 + xn)
            fj3[ny - 2 - i - 1] = (1 - xn) / (1 + xn)

        self.fj3G = cuda.to_device(fj3)
        self.fi3G = cuda.to_device(fi3)
        self.fi2G = cuda.to_device(fi2)
        self.fi1G = cuda.to_device(fi1)
        self.fj2G = cuda.to_device(fj2)
        self.fj1G = cuda.to_device(fj1)
        self.gj2G = cuda.to_device(gj2)
        self.gj3G = cuda.to_device(gj3)
        self.gi2G = cuda.to_device(gi2)
        self.gi3G = cuda.to_device(gi3)

        self.TursG = cuda.to_device(Turs)
        self.LaunchDataG = cuda.to_device(self.LaunchData)
        self.neG = cuda.to_device(self.ne)
        self.AntennaPointsG = cuda.to_device(self.AntennaPoints)
        self.dzG = cuda.to_device(dz)
        self.ezG = cuda.to_device(ez)
        self.hxG = cuda.to_device(hx)
        self.hyG = cuda.to_device(hy)
        self.sxG = cuda.to_device(sx)
        self.sx1G = cuda.to_device(sx1)
        self.sx2G = cuda.to_device(sx2)
        self.omega_peG = cuda.to_device(self.omega_pe)
        self.ihxG = cuda.to_device(ihx)
        self.ihyG = cuda.to_device(ihy)
        threads_per_block = 10
        self.threads_per_block = (threads_per_block,threads_per_block)
        self.blocks_per_grid = (int(nx/threads_per_block+1),int(ny/threads_per_block+1))
        # self.blocks_per_grid = math.ceil(nx*ny / self.threads_per_block)
    def MTANHSi(self,a,x):
        A = a[0]
        B = a[1]
        alpha = a[2]
        x0 = a[3]
        w = a[4]
        y = 0
        z = (x0 - x) / w
        mtanh = ((1 + alpha * z) * np.exp(z) - (1 - 0.04 * z) * np.exp(-z)) / (np.exp(z) + np.exp(-z))
        y = (A * mtanh + B) * 10**19
        return y
    def Caune(self,x,y):
        stpot = 0.25
        ped_pos = 2.25
        h = 4
        w = 0.06
        slope1 = 0.03
        slope2 = -0.1
        a0=[h / 2 ,h / 2 ,slope1 ,ped_pos ,w ,slope2]
        ny = np.size(y)
        nx = np.size(x)
        ne = np.zeros((nx,ny))
        iter1 = np.linspace(0,ny-1,ny).astype(np.int32)
        iter2 = np.linspace(0, nx-1, nx).astype(np.int32)
        for j in iter1:
            if y[j] >= stpot:
                for i in iter2:
                    ne[i,j] = self.MTANHSi(a0, 2.35 - y[j] + stpot)
        return ne


    def omega_pecal(self,ne):
        [rows,cols] = ne.shape
        e = 1.602*10**(-19)
        me = 9.1096*10**(-31)
        epsilon0 = 8.85 * 10**(-12)
        # omega_pe = np.zeros((rows,cols))
        omega_pe = ne*e**2/epsilon0 / me
        return omega_pe
    def setAntenna(self,Ax,waveGL,waveGW,AL,ddx,ddy):
        AntennaPoints = np.zeros((10000,2)).astype(np.int32)
        waveguidthickness = 2
        w = waveGW
        idx1 = int((Ax - w / 1) / ddx) - waveguidthickness - 1
        idx2 = int((Ax - w / 1) / ddx) - 1
        idx3 = int( waveGL/ ddy)
        idx0 = 0
        iters1 = np.linspace(idx1,idx2,idx2-idx1+1).astype(np.int32)
        iters2 = np.linspace(0, idx3-1, idx3).astype(np.int32)
        for i in iters1:
            for j in iters2:
                AntennaPoints[idx0,0] = i
                AntennaPoints[idx0,1] = j
                idx0 += 1
        idx1 = int((Ax + w / 1) / ddx) - 1
        idx2 = int((Ax + w / 1) / ddx) + waveguidthickness - 1
        idx3 = int(waveGL / ddy)
        iters1 = np.linspace(idx1,idx2,idx2-idx1+1).astype(np.int32)
        iters2 = np.linspace(0, idx3-1, idx3).astype(np.int32)
        for i in iters1:
            for j in iters2:
                AntennaPoints[idx0,0] = i
                AntennaPoints[idx0,1] = j
                idx0 += 1
        xm1 = Ax + w
        ym1 = waveGL
        xm2 = Ax - w
        ym2 = waveGL
        ap = 80 * np.pi / 180
        idx1 = int(ym1 / ddy) - 1
        idx2 = int((ym1 + AL) / ddy) - 1
        iters1 = np.linspace(idx1,idx2,idx2-idx1+1).astype(np.int32)
        for j in iters1:
            y1 = (j + 1) * ddy
            x1 = (y1 - ym1) / math.tan(ap) + xm1
            idx3 = int(x1 / ddx) - 1
            idx4 = int(x1 / ddx) + waveguidthickness - 1
            iters2 = np.linspace(idx3, idx4, idx4 - idx3 + 1).astype(np.int32)
            for i in iters2:
                AntennaPoints[idx0,0] = i
                AntennaPoints[idx0,1] = j
                idx0 += 1
        for j in iters1:
            y2 = (j+1) * ddy
            x2 = xm2 - (y2 - ym2) / math.tan(ap)
            idx3 = int(x2 / ddx) - waveguidthickness-1
            idx4 = int(x2 / ddx)-1
            iters2 = np.linspace(idx3, idx4, idx4 - idx3 + 1).astype(np.int32)
            for i in iters2:
                AntennaPoints[idx0,0] = i
                AntennaPoints[idx0,1] = j
                idx0 += 1
        annx = int(idx0)
        AntennaPointss = np.zeros((annx,2)).astype(np.int32)
        AntennaPointss = AntennaPoints[0:annx,:]
        return AntennaPointss

    def CopyDatatoGpu(self):
        pass

    def cauphase(self):
        [rowr,colr] = np.shape(self.ReceiveDatas)
        y1 = np.sum(self.ReceiveDatas,axis=0)/rowr
        self.ReceiveData = y1
        y2 = self.LaunchData
        row = np.size(self.t)
        fs = 1/(self.t[1]-self.t[0])
        T = 1/self.f
        widx = int(T*fs)
        tphase = np.zeros(row).astype(np.float64)
        phase = np.zeros(row).astype(np.float64)
        ReflI = np.zeros(row).astype(np.float64)
        ReflQ = np.zeros(row).astype(np.float64)
        j=0
        i=0
        ysum = y1+y2
        # plt.plot(self.t[0:np.size(y1)],ysum)
        # plt.show()
        y1 = np.abs(y1)
        y2 = np.abs(y2)
        ysum = np.abs(ysum)
        while i < row-widx:
            I1 = math.pow(np.max(y1[i:i+widx]),2)
            # I1 = np.square(np.max(y1[i:i+widx]))
            I2 = math.pow(np.max(y2[i:i + widx]), 2)
            # I2 = np.square(np.max(y2[i:i + widx]))
            Isum = math.pow(np.max(ysum[i:i + widx]), 2)
            # Isum = np.square(np.max(ysum[i:i + widx]))
            if I1 ==0:
                phase[j] = 0
                # I1 = 1*10**(-10)
            else:
                cosph = (Isum - I1 - I2) / (2 * math.pow(I1, 0.5) * math.pow(I2, 0.5))
                if cosph >1:
                    cosph=1
                if cosph < -1:
                    cosph = -1
                phase[j] = math.acos(cosph)
                ReflI[j] = I1 * math.cos(math.acos(cosph))
                ReflQ[j] = I1 * math.sin(math.acos(cosph))
            # cosph = (Isum - I1 - I2) / (2 * np.sqrt(I1)*np.sqrt(I2))
            # cosph = round(cosph, 8)
            # print(y1[i:i+widx])
            # print(I1)
            # print(y2[i:i+widx])
            # print(I2)
            # print(ysum[i:i+widx])
            # print(Isum)
            # print(Isum - I1 - I2)
            # print(2 * np.sqrt(I1)*np.sqrt(I2))
            # print('\n')
            tphase[j] = self.t[i]
            i = i + widx - 1
            j=j+1
        self.tphase = tphase[0:j]
        self.phase = phase[0:j]
        self.ReflI = ReflI[0:j]
        self.ReflQ = ReflQ[0:j]

    def OneDtoTwoD(self,n):
        j = int(n/self.ny)
        i = n % self.ny
        return [i,j]
    def SaveData(self):
        ez = self.ezG.copy_to_host()
        phase = self.phase
        tphase = self.tphase
        LaunchData = self.LaunchData
        ReceiveData = self.ReceiveData
        ReflI = self.ReflI
        ReflQ = self.ReflQ
        t = self.t
        R = self.R
        Z = self.Z
        datas = {
            'ez':ez,
            'phase':phase,
            'tphase':tphase,
            'LaunchData':LaunchData,
            'ReceiveData':ReceiveData,
            'ReflI':ReflI,
            'ReflQ':ReflQ,
            't':t,
            'R':R,
            'Z':Z
        }

        curtrace = os.getcwd()
        if os.path.isdir(curtrace+'/result'):
            sio.savemat('result/FullData.mat',datas)
        else:
            os.makedirs(curtrace+'/result')
            sio.savemat('result/FullData.mat', datas)
    def Gpurun(self):
        print('....../开始迭代/.......')
        start = time.time()
        T = 0
        t0=150
        nx = self.nx
        ny = self.ny
        yatp = 50 * self.ddy
        idx1 = int((self.xatp - self.WGw / 1) / self.ddx)
        idx2 = int((self.xatp + self.WGw / 1) / self.ddx)
        idx3 = int((yatp) / self.ddy)
        idx11 = int((self.Ax - self.WGw / 1) / self.ddx)
        idx22 = int((self.Ax + self.WGw / 1) / self.ddx)
        self.ReceiveDatas = np.zeros((idx22-idx11+1,self.nsteps)).astype(np.float64)
        self.ReceiveDatasG = cuda.to_device(self.ReceiveDatas)

        n_addsource = int(idx2-idx1+1)
        array = np.zeros((nx,ny))
        blockDim_addS = (2, 1)
        gridDim_addS=(int(n_addsource/2+1+70),int(1/1))
        threds_setAnt = 32
        blockDim_setAnt = (threds_setAnt, 1)
        gridDim_setAnt=(int(self.rowAn/threds_setAnt+1),int(1/1))
        Vt = self.Vt
        self.t = np.zeros(self.nsteps)

        for i in range(1,self.nsteps+1):
            # print(i)
            if not i % 1000:
                print('迭代次数：',i)
            T = T + 1
            self.t[i-1] = (i-1) * self.dt
            Dne[self.blocks_per_grid,self.threads_per_block](self.omega_peG,self.neG,Vt, i, self.dt, self.ddx, self.ddy, self.nx, self.ny,self.TursG)
            cuda.synchronize()
            Caudz[self.blocks_per_grid,self.threads_per_block](self.dzG,self.hyG,self.hxG,self.nx, self.ny, self.gi3G, self.gj3G, self.gi2G, self.gj2G, self.ddx, self.ddy)
            cuda.synchronize()
            AddSource[gridDim_addS,blockDim_addS](self.dzG,n_addsource,idx1,idx3,self.ddx,self.xatp,self.f,self.dt,t0,T,self.WGw,self.LaunchDataG)
            cuda.synchronize()
            GetReceiveData[gridDim_addS,blockDim_addS](self.ezG,idx11, idx22, idx3,self.ReceiveDatasG,i)
            cuda.synchronize()

            setAntennaBoundary[gridDim_setAnt,blockDim_setAnt](self.AntennaPointsG,self.rowAn, self.dzG)
            cuda.synchronize()
            Cauez[self.blocks_per_grid, self.threads_per_block](self.ezG, self.dzG, self.sxG, self.sx1G, self.sx2G, self.omega_peG, self.dt, self.vc,self.nx, self.ny)
            cuda.synchronize()
            CauHxy[self.blocks_per_grid, self.threads_per_block](self.ezG, self.ihxG, self.ihyG, self.hxG, self.hyG, self.fj3G, self.fj2G, self.fi1G, self.fi3G, self.fi2G, self.fj1G, self.ddx, self.ddy, self.nx, self.ny)
            cuda.synchronize()
            # plt.figure(1)
            # dz0 = self.dzG.copy_to_host()
            # # dz0 = np.reshape(dz0, (nx, ny))
            # plt.imshow(dz0.T)
            # plt.show()

        end = time.time()
        print('用时：',end-start)
        self.ReceiveDatas = self.ReceiveDatasG.copy_to_host()
        self.LaunchData = self.LaunchDataG.copy_to_host()
        self.cauphase()
        self.SaveData()


        # ez0 = self.ezG.copy_to_host()
        # sio.savemat('ez.mat', {'ez': ez0})
        # # dz0 = np.reshape(dz0,(nx,ny))
        # # xv, yv = np.meshgrid(self.R, self.Z)
        # # x =np.linspace(0,nx*ny,nx*ny)
        # plt.imshow(ez0)
        # # plt.plot(self.tphase,self.phase)
        # plt.show()


f = 30*10**9
c = 3*10**8
nsteps=36000
kx = 200
yc = 0.32
Vt = 1 * c / 2
# Turs = np.linspace(40,50,2).astype(np.float64)
Turs = np.array([50],dtype=np.float64)
# Turs = np.array([50,30,45,40,35],dtype=np.float64)
# Turs = np.array([50,50*(2/3),50*(1/3)],dtype=np.float64)
Test=CauFullwave()
Test._initparameters(f,nsteps,Turs)
Test.Gpurun()
