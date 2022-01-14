#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream> 
#include <iomanip>
#include <MathToolsZz.h>
#include <windows.h>
#include <matlab_plot.h>
#include <opmat.h>
using namespace std;
#define getarrlen1(arr) sizeof(arr)/sizeof(arr[0])
//#define pi 3.1415926535898
//#define pi acos(-1.0)
const double pi = acos(-1.0);
matlab_plot plot;

using namespace std;

void wrtxt1d(double* z, char* fie)
{
	const int nx = _msize(z) / sizeof(*z);
	ofstream outfile(fie);
	if (!outfile) {
		cout << "Unable to open otfile";
		exit(1); // terminate with error 
	}
	for (int i = 0; i < nx; i++) {
		outfile << z[i] << "\t";
	}
	outfile.close();
}
void wrtxt(double** z, char* fie)
{
	const int nx = _msize(z) / sizeof(*z);
	const int ny = _msize(z[0]) / sizeof(*z[0]);
	ofstream outfile(fie);
	//FILE* fp;
	if (!outfile) {
		cout << "Unable to open otfile";
		exit(1); // terminate with error 
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++)
		{
			outfile << z[i][j] << "\t";
			//outfile << setiosflags(ios::fixed) << setprecision(4)<<z[i][j] << "\t";
		}
		outfile << "\n";
	}
	outfile.close();

}
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
double* omega_pe(T& ne)
{
	double e = 1.602 * pow(10, -19);
	double me = 9.1096 * pow(10, -31);
	double epsilon0 = 8.85 * pow(10, -12);
	int N = _msize(ne) / sizeof(*ne);
	double* omega_pe = new double[N];
	double ls;
	for (int i = 0; i < N; i++)
	{
		ls = ne[i] * pow(e, 2) / epsilon0 / me;
		omega_pe[i] = pow(ls, 0.5);
	}
	return omega_pe;
}
template<typename T>
double* omega_ce(T& x, double b0, double x0)
{
	double e = 1.602 * pow(10, -19);
	double me = 9.1096 * pow(10, -31);
	int N = _msize(x) / sizeof(*x);
	double* omega_ce = new double[N];
	double B;
	for (int i = 0; i < N; i++)
	{
		B = b0 * x0 / x[i];
		omega_ce[i] = e * B / me;
	}
	return omega_ce;
}

int main()
{
	plot.open();
	DWORD start, end;
	MathToolsZz matht;
	start = GetTickCount64();
	const double c = 3 * pow(10, 8);
	const double t0 = 150.0;
	const double f = 60 * pow(10, 9);
	const double lambda = c / f;

	const double x_div = 15;
	const double y_div = 3;
	const double ddx = lambda / x_div;//the space seperation, should be 10 times larger than the wave length+
	const double dt = ddx / c / 2;//the time seperation, should be dt = ddx/sqrt(n)/c, here use dt = ddx/2/c;
	const double ddy = lambda / y_div;
	const double spread = lambda / ddx;
	const int bloc = 8;//the location of the wave injection
	const int bwidth = 10;//the width of injected beam
	const int alpha = 100;//the angle between injection beam and the mid-plane
	const int Pulse = 0;
	const int Sinwave = 1;

	const double x_start = 1.9;
	const double x_end = 2.35;
	double* x = matht._linspace(x_start, x_end, ddx);

	const double ped_pos = 2.25;
	const double h = 4;
	const double w = 0.06;
	const double slope1 = 0.03;
	const double slope2 = -0.1;
	double a0[6] = { h / 2 ,h / 2 ,slope1 ,ped_pos ,w ,slope2 };
	double* ne = MTANH(a0, x);
	const int nx = matht.length(x);
	const int ny = 300;
	const int kstart = 1-1;//start location of the plasma material
	const int kend = nx-1;
	const double B0 = 2.0;
	const double x0 = 1.88;//location of the magnetic axis
	double* omega_ce_1D = omega_ce(x, B0, x0);
	double* omega_pe_1D = omega_pe(ne);
	double** omega_pe = matht.CreatMatrix2d(nx, ny, 0);
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			omega_pe[j][i] = omega_pe_1D[j];
		}
	}
	const double nsteps = 3000;
	double** ex = matht.CreatMatrix2d(nx, ny, 0);
	double** ey = matht.CreatMatrix2d(nx, ny, 0);
	double** dx = matht.CreatMatrix2d(nx, ny, 0);
	double** dy = matht.CreatMatrix2d(nx, ny, 0);
	double** hz = matht.CreatMatrix2d(nx, ny, 0);
	//auxiliary parameters
	double** s1 = matht.CreatMatrix2d(nx, ny, 0);
	double** s2 = matht.CreatMatrix2d(nx, ny, 0);
	double** is1 = matht.CreatMatrix2d(nx, ny, 0);
	double** is2 = matht.CreatMatrix2d(nx, ny, 0);
	double** s1m = matht.CreatMatrix2d(nx, ny, 0);
	double** s1mm = matht.CreatMatrix2d(nx, ny, 0);
	double** s2m = matht.CreatMatrix2d(nx, ny, 0);
	double** s2mm = matht.CreatMatrix2d(nx, ny, 0);
	double** idx = matht.CreatMatrix2d(nx, ny, 0);
	double** idy = matht.CreatMatrix2d(nx, ny, 0);
	//Calculat the PML parameters
	double* gi1 = matht.CreatMatrix1d(nx, 0);
	double* gi2 = matht.CreatMatrix1d(nx, 1);
	double* gi3 = matht.CreatMatrix1d(nx, 1);
	double* fi2 = matht.CreatMatrix1d(nx, 1);
	double* fi3 = matht.CreatMatrix1d(nx, 1);
	double* gj1 = matht.CreatMatrix1d(ny, 0);
	double* gj2 = matht.CreatMatrix1d(ny, 1);
	double* gj3 = matht.CreatMatrix1d(ny, 1);
	double* fj2 = matht.CreatMatrix1d(ny, 1);
	double* fj3 = matht.CreatMatrix1d(ny, 1);
	const double vc = 5 * pow(10, 8);
	double* A = matht.CreatMatrix1d(nx, 0);
	double* B = matht.CreatMatrix1d(nx, 0);
	double* C = matht.CreatMatrix1d(nx, 0);

	const double D = exp(-2*vc * dt);
	for (int i = 0; i < nx; i++)
	{
		A[i] = pow(omega_pe_1D[i], 2) / omega_ce_1D[i];
		B[i] = exp(-vc * dt) * sin(omega_ce_1D[i] * dt);
		C[i] = exp(-vc*dt)*cos(omega_ce_1D[i] * dt);
	}

	opmat opt;
	/*
	opmat opt;
	char fe[] = "F:/小黑云同步/学术/Fullwave/Fullwave(C++)/X_mode_2D_boundary_Second/hz.mat";
	char pname1[] = "hz";
	double** lsb=opt.readmat2d(fe,pname1);
	cout << lsb[352][76];
	*/
	//Creat the PML
	int npml = 8;
	for (int i = 1; i <= npml; i++)
	{
		double xnum = npml - i;
		double xd = npml;
		double xxn = xnum / xd;
		double xn = 0.33 * pow(xxn, 3);
		gi1[i - 1] = xn / 2;
		gi1[nx - 2 - i - 1] = xn / 2;
		gi2[i - 1] = 1 / (1 + xn);
		gi2[nx - 2 - i - 1] = 1 / (1 + xn);
		gi3[i - 1] = (1 - xn) / (1 + xn);
		gi3[nx - 2 - i - 1] = (1 - xn) / (1 + xn);
		gj1[i - 1] = xn / 2;
		gj1[ny - 2 - i - 1] = xn / 2;
		gj2[i - 1] = 1 / (1 + xn);
		gj2[ny - 2 - i - 1] = 1 / (1 + xn);
		gj3[i - 1] = (1 - xn) / (1 + xn);
		gj3[ny - 2 - i - 1] = (1 - xn) / (1 + xn);
		xxn = (xnum - 0.5) / xd;
		xn = 0.33 * pow(xxn, 3);
		fi2[i - 1] = 1 / (1 + xn);
		fi2[nx - 2 - i - 1] = 1 / (1 + xn);
		fi3[i - 1] = (1 - xn) / (1 + xn);
		fi3[nx - 2 - i - 1] = (1 - xn) / (1 + xn);
		fj2[i - 1] = 1 / (1 + xn);
		fj2[ny - 2 - i - 1] = 1 / (1 + xn);
		fj3[i - 1] = (1 - xn) / (1 + xn);
		fj3[ny - 2 - i - 1] = (1 - xn) / (1 + xn);
	}

	double** check2 = matht.CreatMatrix2d(nx, ny, 0);
	double* check1 = new double[nx];
	int T = 0;
	for (int n = 1; n <= nsteps; n++)
	{
		cout << n << endl;
		T = T + 1;
		double N_lambda = lambda / ddx;
		double signal = 0;
		double tT = 0;
		for (int is = ny / 2 - bwidth-1; is < ny / 2 + bwidth; is++)
		{
			if (Sinwave == 1)
			{
				tT = t0 - T;
				signal = 10 * sin(-2 * pi * f * dt * (tT));
			}
			hz[nx-bloc-1][is] = hz[nx - bloc-1][is] + signal;
		}
		double curl_hx = 0;
		double curl_hy = 0;
		for (int ii = 2 - 1; ii < nx - 1; ii++)
		{
			for (int jj = 2 - 1; jj < ny - 1; jj++)
			{
				curl_hx = (ddx / ddy) * (hz[ii][jj] - hz[ii][jj - 1]);
				idx[ii][jj] = idx[ii][jj] + curl_hx;
				dx[ii][jj] = gj3[jj] * dx[ii][jj] + gj2[jj] * (0.5 * curl_hx + gi1[ii] * idx[ii][jj]);
				curl_hy = -hz[ii][jj] + hz[ii - 1][jj];
				idy[ii][jj] = idy[ii][jj] + curl_hy;
				dy[ii][jj] = gi3[ii] * dy[ii][jj] + gi2[ii] * (0.5 * curl_hy + gj1[jj] * idy[ii][jj]);
			}
		}
		for (int ii = 2 - 1; ii < nx ; ii++)
		{
			for (int jj = 2 - 1; jj < ny; jj++)
			{
				if (ii >= kstart && ii <= kend)
				{
					ex[ii][jj] = dx[ii][jj] - A[ii] * s1m[ii][jj] - A[ii] * vc * dt * is1[ii][jj] - A[ii] * omega_ce_1D[ii] * dt * is2[ii][jj];
					ey[ii][jj] = dy[ii][jj] - A[ii] * s2m[ii][jj] - A[ii] * vc * dt * is2[ii][jj] + A[ii] * omega_ce_1D[ii] * dt * is1[ii][jj];
					s1[ii][jj] = 2 * C[ii] * s1m[ii][jj] - D * s1mm[ii][jj] + B[ii] * ex[ii][jj] * dt;
					is1[ii][jj] = is1[ii][jj] + s1[ii][jj];
					s2[ii][jj] = 2 * C[ii] * s2m[ii][jj] - D * s2mm[ii][jj] + B[ii] * ey[ii][jj] * dt;
					is2[ii][jj] = is2[ii][jj] + s2[ii][jj];
					s1mm[ii][jj] = s1m[ii][jj];
					s1m[ii][jj] = s1[ii][jj];
					s2mm[ii][jj] = s2m[ii][jj];
					s2m[ii][jj] = s2[ii][jj];
				}
				else
				{
					ex[ii][jj] = dx[ii][jj];
					ey[ii][jj] = dy[ii][jj];
				}
			}
		}
		//calculate the Hz field
		for (int ii = 2 - 1; ii < nx-1; ii++)
		{
			for (int jj = 2 - 1; jj < ny-1; jj++)
			{
				hz[ii][jj] = fi3[ii] * fj3[jj] * hz[ii][jj]+ fi2[ii] * fj2[jj] * 0.5 * ((ddx / ddy) * (ex[ii][jj + 1] - ex[ii][jj])- ey[ii + 1][jj] + ey[ii][jj]);
			}
		}

		plot.plot(x,x);
		//char fie[] = "hz_c++.txt";
		//wrtxt((double**)dy, fie);
		//char fie1[] = "hz_c1++.txt";
		//wrtxt((double**)ey, fie1);
		//char fie2[] = "hz_c2++.txt";
		//wrtxt((double**)hz, fie2);
		//system("pause");
	}
	end = GetTickCount64() - start;
	cout << end << endl;

	char pname[] = "hzc";
	opt.writemat2d(hz, pname);

	char fie[] = "hz_c++.txt";
	//wrtxt((double**)hz, fie);
	system("pause");
	return 0;
}