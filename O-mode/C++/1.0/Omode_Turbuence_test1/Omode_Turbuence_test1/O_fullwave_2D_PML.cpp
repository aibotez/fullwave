#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
//#include "engine.h" //添加MATLAB引擎头文件
//#include "matlab_plot.h"//自定义matlab画图对象
#include <iomanip>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <assert.h>
#include <sstream>
#include <assert.h>
#include <windows.h>
#include <opmat.h>

#include "Efitdensity.h"


using namespace std;
#define getarrlen1(arr) sizeof(arr)/sizeof(arr[0])
//const double pi = M_PI;
const double pi = acos(-1.0);


/*
double** R;
double** Z;
double** ne;
int rows = 0;
int cols = 0;
void read_data()
{
	ifstream infile_R, infile_Z, infile_ne;
	infile_R.open("R.txt");   //将文件流对象与文件连接起来 
	assert(infile_R.is_open());   //若失败,则输出错误消息,并终止程序运行 
	infile_Z.open("Z.txt");
	assert(infile_Z.is_open());
	infile_ne.open("ne.txt");
	assert(infile_ne.is_open());

	string s,sz,sne;
	double b = 0;
	double bz = 0;
	double bne = 0;
	int st = 0;
	int col = 1;
	int row = 0;
	//vector<vector<double>> data;
	vector<vector<double>> data_R;
	vector<vector<double>> data_z;
	vector<vector<double>> data_ne;
	vector<double> vec10;
	vector<double> vec20;
	vector<double> vec30;
	while (getline(infile_R, s))
	{

		getline(infile_Z, sz);
		getline(infile_ne, sne);


		if (col == 1)
		{
			for (int idx = 0; idx < s.length(); idx++)
			{
				if (s[idx] == '\t')
				{
					col = col + 1;
				}
			}
		}
		stringstream stream, stream_z, stream_ne;
		stream << s;
		stream_z << sz;
		stream_ne << sne;
		for (int j = 0; j < col; j++)
		{
			stream >> b;
			vec10.push_back(b);
			stream_z >> bz;
			vec20.push_back(bz);
			stream_ne >> bne;
			vec30.push_back(bne);
		}
		data_R.push_back(vec10);
		vec10.clear();
		data_z.push_back(vec20);
		vec20.clear();
		data_ne.push_back(vec30);
		vec30.clear();

		row = row + 1;
	}
	infile_R.close();
	infile_Z.close();
	infile_ne.close();
	rows = row;
	cols = col;
	R = new double* [row];
	Z = new double* [row];
	ne = new double* [row];
	for (int i = 0; i <row; i++)
	{
		R[i] = new double[col];
		Z[i] = new double[col];
		ne[i] = new double[col];
		for (int j = 0; j < col; j++)
		{
			R[i][j] = data_R[i][j];
			Z[i][j] = data_z[i][j];
			ne[i][j] = pow(10,19)*data_ne[i][j];
		}
		//memcpy(p[i], &data[i], col * sizeof(double));
	}
	//memcpy(p, &data[0], row * col * sizeof(double));

}
*/
double** omega_pecal(double **ne)
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
			ne1 = ne[i][j] * pow(10, 19);
			ls = ne1 * pow(e, 2) / epsilon0 / me;
			omega_pe[i][j] = pow(ls, 0.5);
		}
	}
	return omega_pe;
}

void wrtxt1d(double* z, int nx)
{
	ofstream outfile("ls.txt");
	if (!outfile) {
		cout << "Unable to open otfile";
		exit(1); // terminate with error 
	}
	for (int i = 0; i < nx; i++) {
		outfile << z[i] << "\t";
	}
	outfile.close();
}
void wrtxt(double** z, int nx, int ny)
{
	ofstream outfile("ez_c++.txt");
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

int max_num(double **a,int xl,int yl)
{
	int i, j, max;
	max = a[0][0];
	for (i = 0; i < xl; i++)  
		for (j = 0; j < yl; j++)
			if (a[i][j] > max) max = a[i][j];
	return max;
}
int main()
{
	DWORD start, end;
	double c = 3 * pow(10, 8);
	double f = 30 * pow(10, 9);
	int Pulse = 0;
	int Sinwave = 1;
	char fie[] = "F:/小黑云同步/学术/Fullwave/X_80872_6.117_nh.txt";
	int shot = 80872;
	double time = 6;
	double** R;
	double** Z;
	double** ne;
	double dR;
	double dZ;
	Efitdensity efi;
	efi.efitdensity(fie, shot, time, f,R, Z, ne, dR, dZ);
	
	//read_data();
	const int nx = _msize(R) / sizeof(*R);
	const int ny = _msize(R[0]) / sizeof(*R[0]);
	//const double dR = abs(R[1][0] - R[0][0]);
	//const double dZ = abs(Z[0][1] - Z[0][0]);
	const int kstart = 0;
	const int kend = nx - 100 - 1;
	const double ddx = dR;
	const double ddy = dZ;
	const double dt = ddx / c / 2;
	double t0 = 150;
	const double lambda = c / f;
	const double spread = lambda / ddx;
	int nsteps = 3000;
	double signal =0;
	const double vc = 5 * pow(10, -3);
	const double vt = exp(-vc * dt);
	//double T1 = 1;
	//signal = 5 * sin(-2 * pi * f * dt * (t0 - T1));
	//cout << -2 * pi * f * dt * (t0 - T1)<<" "<<t0-T1<<" "<<dt*f*pi;
	//system("pause");

	double** omega_pe = new double* [nx];
	double** ez = new double* [nx];
	double** dz = new double* [nx];
	double** hx = new double* [nx];
	double** hy = new double* [nx];
	double** sx = new double* [nx];
	double** sx1 = new double* [nx];
	double** sx2 = new double* [nx];
	double** ihx = new double* [nx];
	double** ihy = new double* [nx];


	double* gi2 = new double[nx];
	double* gi3 = new double[nx];
	double* fi1 = new double[nx];
	double* fi2 = new double[nx];
	double* fi3 = new double[nx];
	double* gj2 = new double[ny];
	double* gj3 = new double[ny];
	double* fj1 = new double[ny];
	double* fj2 = new double[ny];
	double* fj3 = new double[ny];
	for (int i = 0; i < nx; i++)
	{
		omega_pe[i] = new double[ny];
		ez[i] = new double[ny];
		dz[i] = new double[ny];
		hx[i] = new double[ny];
		hy[i] = new double[ny];
		sx[i] = new double[ny];
		sx1[i] = new double[ny];
		sx2[i] = new double[ny];
		ihx[i] = new double[ny];
		ihy[i] = new double[ny];
	}
	for (int i = 0; i < nx; i++)
	{
		gi2[i] = 1;
		gi3[i] = 1;
		fi1[i] = 0;
		fi2[i] = 1;
		fi3[i] = 1;
		for (int j = 0; j < ny; j++)
		{
			ez[i][j] = 0;
			dz[i][j] = 0;
			hx[i][j] = 0;
			hy[i][j] = 0;
			sx[i][j] = 0;
			sx1[i][j] = 0;
			sx2[i][j] = 0;
			ihx[i][j] = 0;
			ihy[i][j] = 0;
			if (i == 0)
			{
				gj2[j] = 1;
				gj3[j] = 1;
				fj1[j] = 0;
				fj2[j] = 1;
				fj3[j] = 1;
			}

		}
	}

	omega_pe = omega_pecal((double**)ne);
	//cout << omega_pe[0][0];
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
	//wrtxt1d((double*)fj3, ny);
	//system("pause");

	start = GetTickCount64();
	double T = 0;
	for (int n = 1; n <= nsteps; n++)
	{
		cout << n << endl;
		T = T + 1;
		for (int ii = 1; ii < nx; ii++)
		{
			for (int jj = 1; jj < ny; jj++)
			{
				dz[ii][jj] = gi3[ii] * gj3[jj] * dz[ii][jj] + 0.5 * gi2[ii] * gj2[jj] * (hy[ii][jj] - hy[ii - 1][jj] - (ddx / ddy) * ((hx[ii][jj] - hx[ii][jj - 1])));
			}
		}
		if (Pulse == 1)
		{
			signal = -2.0 * ((t0 - T) / spread) * exp(-1 * pow(((t0 - T) / spread), 2));
		}
		if (Sinwave == 1)
		{
			signal = 5 * sin(-2 * pi * f * dt * (t0 - T));
		}
		for (int is = floor(ny / 2)+100; is <= floor(ny / 2) + 150; is++)
		{
			dz[nx -400 - 1][is - 1] = dz[nx -10 - 1][is - 1] + signal;
			//cout << dz[x_length - 9 - 1][is - 1] << '\t';
		}
		//cout << endl;
		for (int ii = 1; ii < nx; ii++)
		{
			for (int jj = 1; jj < ny; jj++)
			{
				if (ii >= kstart && ii <= kend)
				{
					ez[ii][jj] = dz[ii][jj] - sx[ii][jj];
					sx[ii][jj] = (1 + exp(-vc * dt)) * sx1[ii][jj] - exp(-vc * dt) * sx2[ii][jj] + (pow(omega_pe[ii][jj], 2) * dt / vc) * (1 - exp(-vc * dt)) * ez[ii][jj];
					sx2[ii][jj] = sx1[ii][jj];
					sx1[ii][jj] = sx[ii][jj];
				}
				else
				{
					ez[ii][jj] = dz[ii][jj];
				}
			}
		}
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
	end = GetTickCount64() - start;
	cout << end << endl;
	//wrtxt((double**)ez, nx, ny);
	opmat opt;
	char pname[] = "ezc";
	opt.writemat2d(ez, pname);
	system("pause");
	return 0;
}

