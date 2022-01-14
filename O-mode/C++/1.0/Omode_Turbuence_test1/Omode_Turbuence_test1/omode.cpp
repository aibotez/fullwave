#include<iostream>
#include <windows.h>
#include <math.h>
#include <string>
#include <MathToolsZz.h>
#include <opmat.h>
#include <matlab_plot.h>
#include<cmath>

using namespace std;
#define getarrlen1(arr) sizeof(arr)/sizeof(arr[0])
//const double pi = M_PI;
const double pi = acos(-1.0);
MathToolsZz matht;
matlab_plot plot;
opmat opt;

template<typename T>
double* MTANH(T& a, double* x)
{
	const int n = _msize(x) / sizeof(*x);
	double A = a[0];
	double B = a[1];
	double alpha = a[2];
	double x0 = a[3];
	double w = a[4];
	double* y = new double[n];
	double z;
	double mtanh;
	for (int i = 0; i < n; i++)
	{
		z = (x0 - x[i]) / w;
		mtanh = ((1 + alpha * z) * exp(z) - (1 - (0.04) * z) * exp(-z)) / (exp(z) + exp(-z));
		y[i] = (A * mtanh + B) * pow(10, 19);
	}
	return y;
}
template<typename T>
double MTANHSi(T& a, double x)
{
	double A = a[0];
	double B = a[1];
	double alpha = a[2];
	double x0 = a[3];
	double w = a[4];
	double y = 0;
	double z;
	double mtanh;
	z = (x0 - x) / w;
	mtanh = ((1 + alpha * z) * exp(z) - (1 - (0.04) * z) * exp(-z)) / (exp(z) + exp(-z));
	y = (A * mtanh + B) * pow(10, 19);
	return y;
}
double** Caune(double* x,double* y,double** ne)
{
	double stpot = 0.25;
	const double ped_pos = 2.25;
	const double h = 4;
	const double w = 0.06;
	const double slope1 = 0.03;
	const double slope2 = -0.1;

	double a0[6] = { h / 2 ,h / 2 ,slope1 ,ped_pos ,w ,slope2 };
	int ny = _msize(y) / sizeof(*y);
	int nx = _msize(x) / sizeof(*x);
	int idxstart=0;
	for (int j = 0; j < ny; j++)
	{
		if (y[j] >= stpot)
		{
			for (int i = 0; i < nx; i++)
			{
				ne[i][j] = MTANHSi(a0, 2.35-y[j]+ stpot);
				
			}
		}
	}
	//cout << MTANHSi(a0, 2.35 - 0.3535 + stpot) << endl;;
	return ne;
}
double** omega_pecal(double** ne)
{
	const int rows = _msize(ne) / sizeof(*ne);
	const int cols = _msize(ne[0]) / sizeof(*ne[0]);
	double e = 1.602 * pow(10, -19);
	double me = 9.1096 * pow(10, -31);
	double epsilon0 = 8.85 * pow(10, -12);
	double** omega_pe = new double* [rows];
	double ls;
	for (int i = 0; i < rows; i++)
	{
		omega_pe[i] = new double[cols];
	}

	double ne1 = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			ne1 = ne[i][j];
			ls = 4*pi*ne1 * pow(e, 2) / epsilon0 / me;
			omega_pe[i][j] = ls;
			//omega_pe[i][j] = pow(ls, 0.5);
		}
	}
	return omega_pe;
}

int** CreatMatrix2d(int row, int col, int initvaule)
{
	//std::cout << row << "-----" << col << std::endl;
	int** vms = new int* [row];
	for (int i = 0; i < row; i++)
	{
		vms[i] = new int[col];
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			vms[i][j] = initvaule;
		}
	}
	//cout << vms[0][0] << "  ..."<<endl;
	return vms;
}
int** setAntenna(double Ax, double waveGL, double waveGW, double AL, double ddx, double ddy)
{
	int** AntennaPoints = CreatMatrix2d(10000, 2, 0);
	const int waveguidthickness = 2;
	double w = waveGW;
	int idx1 = int((Ax - w / 1) / ddx) - waveguidthickness-1;
	int idx2 = int((Ax - w / 1) / ddx)-1;
	int idx3 = int(waveGL / ddy);
	//cout << idx1<<'\t' << idx2 <<'\t'<<idx3 << endl;
//system("pause");
	int idx0 = 0;
	for (int i = idx1; i <= idx2; i++)
	{
		for (int j = 0; j < idx3; j++)
		{
			AntennaPoints[idx0][0] = i;
			AntennaPoints[idx0][1] = j;
			idx0 += 1;
		}
	}
	idx1 = int((Ax + w / 1) / ddx)-1;
	idx2 = int((Ax + w / 1) / ddx) + waveguidthickness-1;
	idx3 = int(waveGL / ddy);

	for (int i = idx1; i <= idx2; i++)
	{
		for (int j = 0; j < idx3; j++)
		{
			AntennaPoints[idx0][0] = i;
			AntennaPoints[idx0][1] = j;
			idx0 += 1;
		}
	}
	double xm1 = Ax + w;
	double ym1 = waveGL;
	double xm2 = Ax - w;
	double ym2 = waveGL;
	double ap = 80 * pi / 180;
	idx1 = int(ym1 / ddy)-1;
	idx2 = int((ym1 + AL) / ddy)-1;
	
	//x=int((Ax - w / 1) / ddx) - 1;
	//y= int(waveGL / ddy);
	for (int j = idx1; j <= idx2; j++)
	{
		double y1 = (j+1) * ddy;
		double x1 = (y1 - ym1) / tan(ap) + xm1;
		idx3 = int(x1 / ddx)-1;
		int idx4 = int(x1 / ddx) + waveguidthickness-1;
		for (int i = idx3; i <= idx4; i++)
		{
			AntennaPoints[idx0][0] = i;
			AntennaPoints[idx0][1] = j;
			idx0 += 1;
		}
		
	}
	for (int j = idx1; j <= idx2; j++)
	{
		double y2 = (j+1) * ddy;
		double x2 = xm2 - (y2 - ym2) / tan(ap);
		idx3 = int(x2 / ddx) - waveguidthickness-1;
		int idx4 = int(x2 / ddx)-1;
		for (int i = idx3; i <= idx4; i++)
		{
			AntennaPoints[idx0][0] = i;
			AntennaPoints[idx0][1] = j;
			idx0 += 1;
		}
	}
	int annx = idx0;
	int** AntennaPointss=CreatMatrix2d(annx, 2, 0);
	for (int i = 0; i < annx; i++)
	{
		AntennaPointss[i][0] = AntennaPoints[i][0];
		AntennaPointss[i][1] = AntennaPoints[i][1];
	}
	return AntennaPointss;
}
int** combinematrix(int** x1, int** x2)
{
	const int row1 = _msize(x1) / sizeof(*x1);
	const int row2 = _msize(x2) / sizeof(*x2);
	int** M = CreatMatrix2d(row1 + row2, 2, 0);
	for (int i = 0; i < row1+row2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (i < row1)
			{
				M[i][j] = x1[i][j];
			}
			else
			{
				M[i][j] = x2[i- row1][j];
			}
		}
	}
	return M;
}
void Dne(double** &net,double** & omega_pe,double* x, double* y, double** ne0, double kx, double vt, int ntp, double dt, double ddx, double ddy)
{
	const int row = _msize(x) / sizeof(*x);
	const int col = _msize(y) / sizeof(*y);
	double t = double(ntp-1) * dt;
	double xi = 0;
	double at = 0;
	double Atur = 0;

	double e = 1.602 * pow(10, -19);
	double me = 9.1096 * pow(10, -31);
	double epsilon0 = 8.85 * pow(10, -12);
	double ls;

	double kxi = 0;
	double nesm = 0;
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			xi = i * ddx;
			at = ne0[i][j] * 1 / 100;

			nesm = 0;

			/*
			for (int turn = 1; turn <= 10; turn++)
			{
				kxi = double(turn * kx) / 10;
				Atur = at * sin(kxi * (xi - vt * t));
				nesm = nesm + Atur;
				ls = (ne0[i][j] + nesm) * pow(e, 2) / epsilon0 / me;
				omega_pe[i][j] = ls;
			}
			*/

			xi = i * ddx;
			//at = ne0[i][j] * 1/100;
			Atur = at * sin(kx * (xi - vt * t));
			//xi = i * ddx;
			//at = ne0[i][j] * 1 / 100;
			nesm = nesm + Atur;


			//kxi = 1 * kx / 3;
			//Atur = at * sin(kxi * (xi - vt * t));
			//nesm = nesm + Atur;
			//kxi = 2 * kx / 3;
			//Atur = at * sin(kxi * (xi - vt * t));
			//nesm = nesm + Atur;
			//kxi = 3 * kx / 3;
			//Atur = at * sin(kxi * (xi - vt * t));
			//nesm = nesm + Atur;


			/*kxi = 4 * kx / 10;
			Atur = at * sin(kxi * (xi - vt * t));
			nesm = nesm + Atur;
			kxi = 5 * kx / 10;
			Atur = at * sin(kxi * (xi - vt * t));
			nesm = nesm + Atur;
			kxi = 6 * kx / 10;
			Atur = at * sin(kxi * (xi - vt * t));
			nesm = nesm + Atur;
			kxi = 7 * kx / 10;
			Atur = at * sin(kxi * (xi - vt * t));
			nesm = nesm + Atur;
			kxi = 8 * kx / 10;
			Atur = at * sin(kxi * (xi - vt * t));
			nesm = nesm + Atur;
			kxi = 9 * kx / 10;
			Atur = at * sin(kxi * (xi - vt * t));
			nesm = nesm + Atur;
			kxi = 10 * kx / 10;
			Atur = at * sin(kxi * (xi - vt * t));
			nesm = nesm + Atur;*/
			ls = (ne0[i][j] + nesm) * pow(e, 2) / epsilon0 / me;
			omega_pe[i][j] = ls;


			/*
			Atur = at * sin(kx * (xi - vt * t));
			//net[i][j] = ne0[i][j] + Atur;
			ls = (ne0[i][j] + Atur) * pow(e, 2) / epsilon0 / me;
			omega_pe[i][j] = ls;*/

			//net[i][j] = ls;
		}
	}
}

double GetMaxAbs(double* x,int i1,int i2)
{
	const int row = _msize(x) / sizeof(*x);
	double maxVaule = 0;
	for (int i = i1; i <= i2; i++)
	{
		if (abs(x[i]) > maxVaule)
		{
			maxVaule = abs(x[i]);
		}
	}
	if (maxVaule == 0)
	{
		maxVaule = 1 * pow(10, -5);
	}
	return maxVaule;
}

void cauphase(double* & tphase1,double* & phase1,double* t,double* y1,double* y2,double* ysum,double f, double*& ReflI1, double*& ReflQ1)
{
	const int row = _msize(t) / sizeof(*t);
	double fs = 1 / (t[1] - t[0]);
	double T = 1 / f;
	int widx = int(T*fs);
	//cout << widx << '\t' << endl;
	double* tphase = matht.CreatMatrix1d(row, 0);
	double* phase = matht.CreatMatrix1d(row, 0);
	double* ReflI = matht.CreatMatrix1d(row, 0);
	double* ReflQ = matht.CreatMatrix1d(row, 0);
	int j = 0;
	for (int i = 0; i < row-widx;)
	{
		double I1 = pow(GetMaxAbs(y1,i,i+widx),2);
		double I2 = pow(GetMaxAbs(y2, i, i + widx), 2);
		double Isum = pow(GetMaxAbs(ysum, i, i + widx), 2);
		double cosph = (Isum - I1 - I2) / (2 * pow(I1,0.5)*pow(I2,0.5));
		//cout << acos(cosph) << endl;
		//cout << phase[j]<<'\t' <<j<< endl;
		phase[j] = acos(cosph);
		ReflI[j] = I1 * cos(acos(cosph));
		ReflQ[j] = I1 * sin(acos(cosph));
		tphase[j] = t[i];
		i = i + widx-1;
		j++;
	}
	tphase1 = matht.CreatMatrix1d(j, 0);
	phase1 = matht.CreatMatrix1d(j, 0);
	ReflI1 = matht.CreatMatrix1d(j, 0);
	ReflQ1 = matht.CreatMatrix1d(j, 0);
	for (int i = 0; i < j; i++)
	{
		tphase1[i] = tphase[i];
		phase1[i] = phase[i];
		ReflI1[i] = ReflI[i];
		ReflQ1[i] = ReflQ[i];

	}

}

void savedata(double* R,double* Z,double** ez,double* launchdata,double* receivedata,double* t,double* phaset,double* phase,double* ReflI, double* ReflQ)
{
	char strR[] = "R";
	opt.writemat1d(R,strR);
	char strZ[] = "Z";
	opt.writemat1d(Z,strZ);
	char strez[] = "ez";
	opt.writemat2d(ez, strez);
	char strlaunch[] = "LaunchData";
	opt.writemat1d(launchdata, strlaunch);
	char strreceive[] = "ReceiveData";
	opt.writemat1d(receivedata, strreceive);
	char strt[] = "t";
	opt.writemat1d(t, strt);
	char strtphase[] = "tphase";
	opt.writemat1d(phaset, strtphase);
	char strphase[] = "phase";
	opt.writemat1d(phase, strphase);
	char strReflI[] = "ReflI";
	opt.writemat1d(ReflI, strReflI);
	char strReflQ[] = "ReflQ";
	opt.writemat1d(ReflQ, strReflQ);
}
int main()
{
	
	DWORD start, end;
	const double c = 3 * pow(10, 8);
	const double f = 33 * pow(10, 9);
	int Pulse = 0;
	int Sinwave = 1;
	double* R;
	double* Z;
	double** ne;
	double dR;
	double dZ;
	double lambda = c / f;
	double k0 = 2 * pi / lambda;
	double ddx = lambda / 20;
	double ddy = lambda / 20;
	double dt = ddx / c / 2;
	R = matht._linspace(0, 0.4, ddx);
	Z = matht._linspace(0, 0.5, ddx);
	const int nx = matht.length(R);
	const int ny = matht.length(Z);

	const int kstart = 0;
	const int kend = nx - 100 - 1;
	double t0 = 150;
	//int nsteps = 3000;// 12000 * 3;//12000*3;//480000
	double signal = 0;
	const double vc = 5 * pow(10, -3);
	const double vt = exp(-vc * dt);
	const double w = 0.01 / 2;
	cout << nx <<'\t'<< ny << endl;
	double** ne0 = matht.CreatMatrix2d(nx, ny, 0);
	ne = Caune(R, Z, ne0);
	//char pname[] = "ne";
	//opt.writemat2d(ne, pname);
	//plot.open();
	//plot.imagesc(R,Z,ne);
	//system("pause");
	double xatp = 0.2;
	double Ax = 0.245;
	double AtL = 0.05;
	int** LaunchAntennaPoint = setAntenna(xatp, 1.5 * AtL, w, 1.3 * AtL, ddx, ddy);
	int** ReceiveAntennaPoint = setAntenna(Ax, 1.5 * AtL, w, 1.3 * AtL, ddx, ddy);
	//int** AntennaPoint = LaunchAntennaPoint;
	int** AntennaPoint = combinematrix(LaunchAntennaPoint, ReceiveAntennaPoint);
	const int rowAn = _msize(AntennaPoint) / sizeof(*AntennaPoint);
	double** omega_pe = omega_pecal(ne);
	double** ez = matht.CreatMatrix2d(nx, ny, 0);
	double** dz = matht.CreatMatrix2d(nx, ny, 0);
	double** hx = matht.CreatMatrix2d(nx, ny, 0);
	double** hy = matht.CreatMatrix2d(nx, ny, 0);
	double** sx = matht.CreatMatrix2d(nx, ny, 0);
	double** sx1 = matht.CreatMatrix2d(nx, ny, 0);
	double** sx2 = matht.CreatMatrix2d(nx, ny, 0);
	double** ihx = matht.CreatMatrix2d(nx, ny, 0);
	double** ihy = matht.CreatMatrix2d(nx, ny, 0);
	//Calculat the PML parameters
	double* gi2 = matht.CreatMatrix1d(nx, 1);
	double* gi3 = matht.CreatMatrix1d(nx, 1);
	double* fi1 = matht.CreatMatrix1d(nx, 0);
	double* fi2 = matht.CreatMatrix1d(nx, 1);
	double* fi3 = matht.CreatMatrix1d(nx, 1);

	double* gj2 = matht.CreatMatrix1d(ny, 1);
	double* gj3 = matht.CreatMatrix1d(ny, 1);
	double* fj1 = matht.CreatMatrix1d(ny, 0);
	double* fj2 = matht.CreatMatrix1d(ny, 1);
	double* fj3 = matht.CreatMatrix1d(ny, 1);

	const double coeff_pml = 0.33;
	int npml = 8;
	for (int i = 1; i <= npml; i++)
	{
		double xnum = npml - i;
		double xd = npml;
		double xxn = xnum / xd;
		double xn = coeff_pml * pow(xxn, 3);
		gi2[i - 1] = 1 / (1 + xn);
		gi2[nx - 1 - i - 1] = 1 / (1 + xn);
		gi3[i - 1] = (1 - xn) / (1 + xn);
		gi3[nx - 1 - i - 1] = (1 - xn) / (1 + xn);
		gj2[i - 1] = 1 / (1 + xn);
		gj2[ny - 1 - i - 1] = 1 / (1 + xn);
		gj3[i - 1] = (1 - xn) / (1 + xn);
		gj3[ny - 1 - i - 1] = (1 - xn) / (1 + xn);
		xxn = (xnum - 0.5) / xd;
		xn = coeff_pml * pow(xxn, 3);
		fi1[i - 1] = xn / 2;
		fi1[nx - 2 - i - 1] = xn / 2;
		fi2[i - 1] = 1 / (1 + xn);
		fi2[nx - 2 - i - 1] = 1 / (1 + xn);
		fi3[i - 1] = (1 - xn) / (1 + xn);
		fi2[nx - 2 - i - 1] = (1 - xn) / (1 + xn);
		fj1[i - 1] = xn / 2;
		fj1[ny - 2 - i - 1] = xn / 2;
		fj2[i - 1] = 1 / (1 + xn);
		fj2[ny - 2 - i - 1] = 1 / (1 + xn);
		fj3[i - 1] = (1 - xn) / (1 + xn);
		fj3[ny - 2 - i - 1] = (1 - xn) / (1 + xn);

	}

	start = GetTickCount64();
	double T = 0;
	double Vt = 1 * c / 20;
	double kx = 50;
	double ftur = kx * Vt / (2 * pi);
	double Ttur = 1 / (Vt * kx);
	double** net = matht.CreatMatrix2d(nx, ny, 0);
	int nsteps = 30000;// 12000 * 3;//12000*3;//480000
	double* t = matht.CreatMatrix1d(nsteps,0);
	double sigtmp = 0;
	double UA = 5;
	double Uofft = 100;
	double Ui = 0;
	double Recesigtmp = 0;
	double* LaunchData = matht.CreatMatrix1d(nsteps, 0);
	double* ReceiveData = matht.CreatMatrix1d(nsteps, 0);
	double* SumData = matht.CreatMatrix1d(nsteps, 0);


	double kxs[] = { 110,40,30,50,50,50,50,50,50 };
	double e = 1.602 * pow(10, -19);
	double me = 9.1096 * pow(10, -31);
	double epsilon0 = 8.85 * pow(10, -12);
	double ls;
	for (int n = 1; n <= nsteps; n++)
	{
		//cout << n << endl;
		if (n % 1000 == 0)
		{
			cout << n << endl;
		}

		T = T + 1;
		t[n-1] = double(n-1) * dt;

		//Dne(net, omega_pe,R, Z, ne0, kx, Vt, n, dt, ddx, ddy);
		

	
		
#pragma omp parallel for
		for (int ii = 0; ii < nx; ii++)
		{
			double ls;
			double kxi = 0;
			double nesm = 0;
			double at = 0;
			double xi = ii * ddx;
			double Atur = 0;
			for (int jj = 0; jj < ny; jj++)
			{
				nesm = 0;

				for (int m = 0; m < 1; m++)
				{
					at = ne[ii][jj] * 1 / 100;
					kxi = kxs[m];
					Atur = at * sin(kxi * (xi - Vt * (n-1) * dt));
					nesm = nesm + Atur;
				}
				ls = 4*pi*(ne[ii][jj] + nesm) * pow(e, 2) / epsilon0 / me;
				omega_pe[ii][jj] = ls;

			}
		}
		
		

		//#pragma omp parallel for
		for (int ii = 1; ii < nx; ii++)
		{
			for (int jj = 1; jj < ny; jj++)
			{
				dz[ii][jj] = gi3[ii] * gj3[jj] * dz[ii][jj] + 0.5 * gi2[ii] * gj2[jj] * (hy[ii][jj] - hy[ii - 1][jj] - (ddx / ddy) * ((hx[ii][jj] - hx[ii][jj - 1])));
			}
		}

	//add source
		double yatp = 50 * ddy;
		int idx1 = int((xatp - w / 1) / ddx);
		int idx2 = (xatp + w / 1) / ddx;
		int idx3 = int((yatp) / ddy);
		sigtmp = 0;
		Recesigtmp = 0;
		double ns = n - 1;
		for (int i = idx1; i <= idx2; i++)
		{
			double pl = i * ddx - xatp;
	
			 Ui= UA * (exp(ns / Uofft) - exp(-ns / Uofft)) / (exp(ns / Uofft) + exp(-ns / Uofft));

			signal = Ui * exp(-pow(pl,2)/ pow(w,2)) * cos(-2 * pi * f * dt * (t0 - T));
			sigtmp = sigtmp + signal;
			//Recesigtmp = Recesigtmp+ez[]
			dz[i-1][idx3-1] = dz[i-1][idx3-1] + signal;
		}
		double sumtmp = idx2 - idx1 + 1;
		LaunchData[n-1] = sigtmp / sumtmp;

		idx1 = int((Ax - w / 1) / ddx);
		idx2 = (Ax + w / 1) / ddx;
		for (int i = idx1; i <= idx2; i++)
		{
			Recesigtmp = Recesigtmp + ez[i - 1][idx3 - 1];
		}
		sumtmp = idx2 - idx1 + 1;
		ReceiveData[n-1] = Recesigtmp / sumtmp;
		SumData[n - 1] = LaunchData[n - 1] + ReceiveData[n - 1];

		for (int i = 0; i < rowAn; i++)
		{
			dz[AntennaPoint[i][0]][AntennaPoint[i][1]] = 0;
		}

		#pragma omp parallel for
		for (int ii = 1; ii < nx; ii++)
		{
			for (int jj = 1; jj < ny; jj++)
			{
				ez[ii][jj] = dz[ii][jj] - sx[ii][jj];
				sx[ii][jj] = (1 + exp(-vc * dt)) * sx1[ii][jj] - exp(-vc * dt) * sx2[ii][jj] + (omega_pe[ii][jj] * dt / vc) * (1 - exp(-vc * dt)) * ez[ii][jj];
				sx2[ii][jj] = sx1[ii][jj];
				sx1[ii][jj] = sx[ii][jj];
				/*
				if (ii >= kstart && ii <= kend)
				{
					ez[ii][jj] = dz[ii][jj] - sx[ii][jj];
					sx[ii][jj] = (1 + exp(-vc * dt)) * sx1[ii][jj] - exp(-vc * dt) * sx2[ii][jj] + (omega_pe[ii][jj] * dt / vc) * (1 - exp(-vc * dt)) * ez[ii][jj];
					sx2[ii][jj] = sx1[ii][jj];
					sx1[ii][jj] = sx[ii][jj];
				}
				else
				{
					ez[ii][jj] = dz[ii][jj];
				}*/
			}
		}
		#pragma omp parallel for
		//calculate the Hy field
		for (int ii = 0; ii < nx - 1; ii++)
		{
			for (int jj = 0; jj < ny - 1; jj++)
			{
				double curl_e = ez[ii][jj] - ez[ii][jj + 1];
				ihx[ii][jj] = ihx[ii][jj] + curl_e;
				hx[ii][jj] = fj3[jj] * (ddx / ddy) * hx[ii][jj] + fj2[jj] * 0.5 * (ddx / ddy) * (curl_e + fi1[ii] * ihx[ii][jj]);
			}
		}
		#pragma omp parallel for
		//calculate the Hx field
		for (int ii = 0; ii < nx - 1; ii++)
		{
			for (int jj = 0; jj < ny - 1; jj++)
			{
				double curl_e = ez[ii][jj] - ez[ii + 1][jj];
				ihy[ii][jj] = ihy[ii][jj] + curl_e;
				hy[ii][jj] = fi3[ii] * hy[ii][jj] - fi2[ii] * (0.5 * curl_e + fj1[jj] * ihy[ii][jj]);
			}
		}


	}

	double* tphase;
	double* phase;
	double* ReflI;
	double* ReflQ;
	cauphase(tphase,phase,t, ReceiveData, LaunchData, SumData,f,ReflI,ReflQ);


	end = GetTickCount64() - start;
	cout << end/1000 << endl;

	double** Ant = matht.CreatMatrix2d(rowAn, 2, 0);
	for (int i = 0; i < rowAn; i++)
	{
		Ant[i][0] = AntennaPoint[i][0];
		Ant[i][1] = AntennaPoint[i][1];
	}
	savedata(R,Z,ez,LaunchData,ReceiveData,t,tphase,phase, ReflI, ReflQ);
	//char pname[] = "AntennaPoint";
	//opt.writemat2d(Ant, pname);
	/*
	char pname[] = "phase";
	opt.writemat1d(phase, pname);
	char pnametphase[] = "tphase";
	opt.writemat1d(tphase, pnametphase);
	char pname1[] = "Ez";
	opt.writemat2d(ez, pname1);
	char pnamet[] = "t";
	opt.writemat1d(t, pnamet);
	char pnamerece[] = "ReceiveData";
	opt.writemat1d(ReceiveData, pnamerece);
	char pnamelau[] = "LaunchData";
	opt.writemat1d(LaunchData, pnamelau);*/

	//ReceiveData
	return 0;
}