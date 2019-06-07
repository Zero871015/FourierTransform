#include "FT.h"
#include <complex>

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
		pFreq[newcnt] = new double[N]; // �ť߸��W�v�}�C
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
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
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
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
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
		pFreq[i] = new double[N]; // �ť߸��W�v�}�C
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
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
			//�s�U�ϳť߸���ƻP��Ƴ���
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
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}


/*-------------------------------------------------------------------------
   Calculate the closest but lower power of two of a number
   twopm = 2**m <= n
   Return TRUE if 2**m == n
*/
void Powerof2(int n, int *m, int *twopm)
{
	if (n <= 1) {
		*m = 0;
		*twopm = 1;
	}

	*m = 1;
	*twopm = 2;
	do {
		(*m)++;
		(*twopm) *= 2;
	} while (2 * (*twopm) <= n);
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

	 Formula: forward
				  N-1
				  ---
			  1   \          - j k 2 pi n / N
	  X(n) = ---   >   x(k) e                    = forward transform
			  N   /                                n=0..N-1
				  ---
				  k=0

	  Formula: reverse
				  N-1
				  ---
				  \          j k 2 pi n / N
	  X(n) =       >   x(k) e                    = forward transform
				  /                                n=0..N-1
				  ---
				  k=0
*/
void FFT(int dir, long m, std::complex <double> x[])
{
	long i, i1, i2, j, k, l, l1, l2, n;
	std::complex <double> tx, t1, u, c;

	/*Calculate the number of points */
	n = 1;
	for (i = 0; i < m; i++)
		n <<= 1;

	/* Do the bit reversal */
	i2 = n >> 1;
	j = 0;

	for (i = 0; i < n - 1; i++)
	{
		if (i < j)
			swap(x[i], x[j]);

		k = i2;

		while (k <= j)
		{
			j -= k;
			k >>= 1;
		}

		j += k;
	}

	/* Compute the FFT */
	c.real(-1.0);
	c.imag(0.0);
	l2 = 1;
	for (l = 0; l < m; l++)
	{
		l1 = l2;
		l2 <<= 1;
		u.real(1.0);
		u.imag(0.0);

		for (j = 0; j < l1; j++)
		{
			for (i = j; i < n; i += l2)
			{
				i1 = i + l1;
				t1 = u * x[i1];
				x[i1] = x[i] - t1;
				x[i] += t1;
			}

			u = u * c;
		}

		c.imag(sqrt((1.0 - c.real()) / 2.0));
		if (dir == 1)
			c.imag(-c.imag());
		c.real(sqrt((1.0 + c.real()) / 2.0));
	}

	/* Scaling for forward transform */
	if (dir == 1)
	{
		for (i = 0; i < n; i++)
			x[i] /= n;
	}
	return;
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
	  the dimensions are not powers of 2
*/
void FFT2D(std::complex<double> **c, int nx, int ny, int dir)
{
	int i, j;
	int m, twopm;

	Powerof2(nx, &m, &twopm);

	/* Transform the rows */
	std::complex <double> *x = new std::complex <double>[nx];

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			x[i].real(c[i][j].real());
			x[i].imag(c[i][j].imag());
		}
		FFT(dir, m, x);
		for (i = 0; i < nx; i++) {
			c[i][j].real(x[i].real());
			c[i][j].imag(x[i].imag());
		}
	}
	delete[] x;

	/* Transform the columns */
	x = new std::complex <double>[nx];

	Powerof2(ny, &m, &twopm);
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			x[j].real(c[i][j].real());
			x[j].imag(c[i][j].imag());
		}
		FFT(dir, m, x);
		for (j = 0; j < ny; j++) {
			c[i][j].real(x[j].real());
			c[i][j].imag(x[j].imag());
		}
	}
	delete[] x;
}



void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // �ť߸��W�v�}�C
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	
	std::complex<double> ** input_seq;
	input_seq = new std::complex<double>*[h];
	for (unsigned int j = 0; j < h; j++) {
		input_seq[j] = new std::complex<double>[h];
	}
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; ++j) {
			input_seq[i][j] = std::complex<double>(InputImage[i][j], 0);
		}
	}

	FFT2D(input_seq, w, h, -1);

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(input_seq[i][j].real(), (double) 2.0) + pow(input_seq[i][j].imag(), (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
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
	// M = N �����O��}
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

			// �Ck�Ӱ��@��FFT
			for (int j = 0; j < N; j += k)
			{
				// �ek/2�ӻP��k/2���T����ƭȫ�n��١A
				// �]������٪��@�_���C
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

			// �Ck�Ӱ��@��FFT
			for (int j = 0; j < N; j += k)
			{
				// �ek/2�ӻP��k/2���T����ƭȫ�n��١A
				// �]������٪��@�_���C
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
		pFreq[newcnt] = new double[N]; // �ť߸��W�v�}�C
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------

	std::complex<double> ** input_seq;
	input_seq = new std::complex<double>*[h];
	for (unsigned int j = 0; j < h; j++) {
		input_seq[j] = new std::complex<double>[h];
	}
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; ++j) {
			input_seq[i][j] = std::complex<double>(InputImage[i][j], 0);
		}
	}

	FFT2D(input_seq, w, h, 1);

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(input_seq[i][j].real(), (double) 2.0) + pow(input_seq[i][j].imag(), (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
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
