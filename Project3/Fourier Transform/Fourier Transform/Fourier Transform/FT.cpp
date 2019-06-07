#include "FT.h"
#include "complex.h"

FT::FT()
{
	
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt<M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i<M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j<N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage,M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i<M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag,FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
			//存下反傅立葉實數與虛數部分
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

int rev_bits(unsigned int index, int size)
{
	int rev = 0;
	for (; size > 1; size >>= 1)
	{
		rev = (rev << 1) | (index & 1);
		index >>= 1;
	}
	return rev;
}

Complex* fft1(Complex *data, unsigned int size, int log2n, bool inverse)
{
	double angle, wtmp, wpr, wpi, wr, wi;
	int n = 1, n2;
	double pi2 = 3.14159 * 2.0;
	double scale = 1.0 / size;
	Complex tc;

	for (int k = 0; k < log2n; ++k)
	{
		n2 = n;
		n <<= 1;
		angle = (inverse) ? pi2 / n : -pi2 / n;
		wtmp = sin(0.5*angle);
		wpr = -2.0*wtmp*wtmp;
		wpi = sin(angle);
		wr = 1.0;
		wi = 0.0;

		for (int m = 0; m < n2; ++m) {
			for (unsigned int i = m; i < size; i += n) {
				int j = i + n2;
				tc.real = wr * data[j].real - wi * data[j].imag;
				tc.imag = wr * data[j].imag + wi * data[j].real;
				data[j] = data[i] - tc;
				data[i] += tc;
			}
			wr = (wtmp = wr)*wpr - wi * wpi + wr;
			wi = wi * wpr + wtmp * wpi + wi;
		}
	}
	if (inverse) {
		for (int i = 0; i < n; i++) {
			data[i] *= scale;
		}
	}
	return data;
}

Complex** fft2(Complex **data, int r, int c, bool inverse)
{
	Complex *row = new Complex[r];
	Complex *column = new Complex[c];
	int log2w = log2(r);
	int log2h = log2(c);

	// Perform FFT on each row
	for (int y = 0; y < c; ++y) {
		for (int x = 0; x < r; ++x) {
			int rb = rev_bits(x, r);
			row[rb] = data[x][y];
		}
		row = fft1(row, r, log2w, inverse);
		for (int x = 0; x < r; ++x) {
			data[x][y] = row[x];
		}
	}

	// Perform FFT on each column
	for (int x = 0; x < r; ++x) {
		for (int y = 0; y < c; ++y) {
			int rb = rev_bits(y, c);
			column[rb] = data[x][y];
		}
		column = fft1(column, c, log2h, inverse);
		for (int y = 0; y < c; ++y) {
			data[x][y] = column[y];
		}
	}
	return data;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	/*
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			FFT(FreqReal, FreqImag, InputImage, M, N, j, i);
		}
	}
	*/
	Complex ** input_seq;
	input_seq = new Complex*[h];
	for (unsigned int j = 0; j < h; j++) {
		input_seq[j] = new Complex[h];
	}
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; ++j) {
			input_seq[i][j] = Complex(InputImage[i][j], 0);
		}
	}

	fft2(input_seq, h, w, false);

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(input_seq[i][j].real, (double) 2.0) + pow(input_seq[i][j].imag, (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

/*
void FT::FFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;
	
	for (int q = 0; q < N; q++)
	{
		std::complex<double> *x = new std::complex < double >[N] ;
		for (int i = 0; i < M; i++)
		{
			x[i] = InputImage[i][q];
		}

		for (int i = 1, j = 0; i < N; ++i)
		{
			for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
			//      for (int k=N>>1; k>(j^=k); k>>=1) ;
			if (i > j) swap(x[i], x[j]);
			//      if (i<j) swap(x[i], x[j]);
		}



		for (int k = 2; k <= N; k <<= 1)
		{
			float theta = -2.0 * 3.14159 / k;
			std::complex<double> delta_w(cos(theta), sin(theta));

			// 每k個做一次FFT
			for (int j = 0; j < N; j += k)
			{
				// 前k/2個與後k/2的三角函數值恰好對稱，
				// 因此兩兩對稱的一起做。
				std::complex < double > w(1, 0);
				for (int i = j; i < j + k / 2; i++)
				{
					std::complex < double > a = x[i];
					std::complex<double> b = x[i + k / 2] * w;
					x[i] = a + b;
					x[i + k / 2] = a - b;
					w *= delta_w;
				}
			}
		}

		for (int i = 0; i < M; i++)
		{
			if (i > M / 2)
			{
				pFreqReal[u][v] += x[i - M / 2].real();
				pFreqImag[u][v] -= x[i - M / 2].imag();
			}
			else
			{
				pFreqReal[u][v] += x[i + M / 2].real();
				pFreqImag[u][v] -= x[i + M / 2].imag();
			}
		}
		delete[] x;
	}


	for (int q = 0; q < M; q++)
	{
		std::complex<double> *x = new std::complex < double >[M];
		for (int i = 0; i < N; i++)
		{
			x[i] = InputImage[q][i];
		}

		for (int i = 1, j = 0; i < M; ++i)
		{
			for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
			//      for (int k=N>>1; k>(j^=k); k>>=1) ;
			if (i > j) swap(x[i], x[j]);
			//      if (i<j) swap(x[i], x[j]);
		}

		for (int k = 2; k <= N; k <<= 1)
		{
			float theta = -2.0 * 3.14159 / k;
			std::complex<double> delta_w(cos(theta), sin(theta));

			// 每k個做一次FFT
			for (int j = 0; j < N; j += k)
			{
				// 前k/2個與後k/2的三角函數值恰好對稱，
				// 因此兩兩對稱的一起做。
				std::complex < double > w(1, 0);
				for (int i = j; i < j + k / 2; i++)
				{
					std::complex < double > a = x[i];
					std::complex<double> b = x[i + k / 2] * w;
					x[i] = a + b;
					x[i + k / 2] = a - b;
					w *= delta_w;
				}
			}
		}

		for (int i = 0; i < N; i++)
		{
			if (i > M / 2)
			{
				pFreqReal[u][v] += x[i - M / 2].real();
				pFreqImag[u][v] -= x[i - M / 2].imag();
			}
			else
			{
				pFreqReal[u][v] += x[i + M / 2].real();
				pFreqImag[u][v] -= x[i + M / 2].imag();
			}
		}
		delete[] x;
	}
	
	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}
*/

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	/*
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			FFT(FreqReal, FreqImag, InputImage, M, N, j, i);
		}
	}
	*/
	Complex ** input_seq;
	input_seq = new Complex*[h];
	for (unsigned int j = 0; j < h; j++) {
		input_seq[j] = new Complex[h];
	}
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; ++j) {
			input_seq[i][j] = Complex(InputImage[i][j], 0);
		}
	}

	fft2(input_seq, h, w, true);


	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(input_seq[i][j].real, (double) 2.0) + pow(input_seq[i][j].imag, (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
}


void FT::LowpassFilter(double** Real, double** Img, double** filter)
{
}

void FT::HighpassFilter(double** Real, double** Img, double** filter)
{
}
